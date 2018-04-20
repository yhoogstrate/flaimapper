#!/usr/bin/env python

"""FlaiMapper: computational annotation of small ncRNA derived fragments using RNA-seq high throughput data

 Here we present Fragment Location Annotation Identification mapper
 (FlaiMapper), a method that extracts and annotates the locations of
 sncRNA-derived RNAs (sncdRNAs). These sncdRNAs are often detected in
 sequencing data and observed as fragments of their  precursor sncRNA.
 Using small RNA-seq read alignments, FlaiMapper is able to annotate
 fragments primarily by peak-detection on the start and  end position
 densities followed by filtering and a reconstruction processes.
 Copyright (C) 2011-2014:
 - Youri Hoogstrate
 - Elena S. Martens-Uzunova
 - Guido Jenster


 [License: GPL3]

 This file is part of flaimapper.

 flaimapper is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 flaimapper is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 Documentation as defined by:
 <http://epydoc.sourceforge.net/manual-fields.html#fields-synonyms>
"""

import flaimapper
import operator
import logging

from flaimapper.BAMParser import BAMParser
from flaimapper.ncRNAFragment import ncRNAFragment
from flaimapper.utils import sort_frequency_dict
from flaimapper.utils import py2_round


logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)


class MaskedRegion:
    """A masked region is a region masked in the reference genome to
    indicate where ncRNAs are located.
    """
    def __init__(self, region, settings):
        logging.debug("Masked region: " + region[0] + ":" + str(region[1]) + "-" + str(region[2]))

        self.region = region
        self.settings = settings

    def get_median_of_map(self, value_map):
        """
        input:
        {
            21: 2
            23: 6
            24: 3
        }

        output:
            median of [21,21,23,23,23,23,23,23,24,24,24]

        ----

        trick:
        iter1:

        {
            21: 2  -- look at count of lowest key -> 2
            23: 6
            24: 3  -- look at count of larges key -> 3
        }

        remove the key with the lowest count and substract that number from the highest count:
        {   23: 6
            24: 1  }

        iter2:

        {
            23: 6  -- key: 6
            24: 1  -- key: 1
         }

        remove lowest key:
        {
            23: 5
         }

         median = 23


        ******

        in case of equal keys and n>2, remove both:
      { 21:5
        23:4
        24:5 }
        --> { 23:4 }

        ******
        in case of equal keys and n==2, return get_median([k1,k2])

      { 22:8
        26:8 }
        -> median([22,26]) == return (float(k1 + k2)) / 2
        """

        keys = sorted(value_map.keys())
        while len(keys) > 1:
            key_l = keys[0]
            key_t = keys[-1]

            if value_map[key_l] < value_map[key_t]:
                value_map[key_t] -= value_map[key_l]
                del(value_map[key_l])
                keys.remove(key_l)
            elif value_map[key_t] < value_map[key_l]:
                value_map[key_l] -= value_map[key_t]
                del(value_map[key_t])
                keys.remove(key_t)
            elif value_map[key_t] == value_map[key_l]:
                if len(keys) == 2:
                    return (float(keys[0] + keys[1])) / 2
                else:
                    del(value_map[key_l], value_map[key_t])
                    keys.remove(key_l)
                    keys.remove(key_t)

        if len(keys) > 0:
            return keys[0]
        else:
            return None

    def predict_fragments(self):
        def step_01__parse_stats():
            logging.debug("Acquiring statistics")
            n = self.region[2] - self.region[1] + 1  # both zero based; 0-0=0 while that should be 1, so 0-0+1=1

            self_start_positions = [0] * n
            self_stop_positions = [0] * n

            tmp_start_avg_lengths = [{} for x in range(n)]  # [{}] * n makes references instead of copies
            tmp_stop_avg_lengths = [{} for x in range(n)]  # [{}] * n makes references instead of copies

            for read in BAMParser(self.region, self.settings.alignment_file):
                pos_start = read[0] - self.region[1]
                pos_stop = read[1] - self.region[1]

                if pos_start >= 0 and pos_stop >= 0 and pos_start < n and pos_stop < n:
                    len_start = read[1] - read[0]
                    len_stop = read[0] - read[1]

                    self_start_positions[read[0] - self.region[1]] += 1
                    self_stop_positions[read[1] - self.region[1]] += 1

                    if len_start not in tmp_start_avg_lengths[pos_start]:
                        tmp_start_avg_lengths[pos_start][len_start] = 0

                    if len_stop not in tmp_stop_avg_lengths[pos_stop]:
                        tmp_stop_avg_lengths[pos_stop][len_stop] = 0

                    tmp_start_avg_lengths[pos_start][len_start] += 1
                    tmp_stop_avg_lengths[pos_stop][len_stop] += 1

                else:
                    logging.error("Alignment out of bound: (%i,%i) %s:%i-%i" % (pos_start, pos_stop, self.region[0], self.region[1], self.region[2]))

            # Calc medians
            self_start_avg_lengths = []
            self_stop_avg_lengths = []

            for i in range(len(tmp_stop_avg_lengths)):
                avgLenF = self.get_median_of_map(tmp_start_avg_lengths[i])
                avgLenR = self.get_median_of_map(tmp_stop_avg_lengths[i])

                if(avgLenF):
                    avgLenF = int(py2_round(avgLenF + 1))
                if(avgLenR):
                    avgLenR = int(py2_round(avgLenR - 0.5))							# Why -0.5 -> because of rounding a negative number

                self_start_avg_lengths.append(avgLenF)
                self_stop_avg_lengths.append(avgLenR)

            return (self_start_positions,
                    self_stop_positions,
                    self_start_avg_lengths,
                    self_stop_avg_lengths)

        def step_02__find_peaks(plist, drop_cutoff=0.1):
            # Define variables:
            peaks = {}

            previous = 0
            highest = 0
            highestPos = -1

            # Walk over list of [start/stop]-position counts:
            for pos in range(len(plist)):
                current = plist[pos]
                if current > previous:  # and (current > (noise_type_alpha_cutoff/100.0*max(plist)))):
                    if current > highest:
                        highest = current
                        highestPos = pos
                elif current < previous:
                    # if (current < (drop_cutoff*highest)) and (highestPos != -1):
                    # if (current < (100.0*drop_cutoff*highest)) and (highestPos != -1):
                    # if (current < (10.0*highest)) and (highestPos != -1):
                    if (drop_cutoff * current < highest) and (highestPos != -1):
                        peaks[highestPos] = highest
                        # highestPos = -1
                        highest = 0

                previous = current

            return peaks

        def step_03__smooth_filter_peaks(plist):
            """Smooth filtering
            """

            psorted = sorted(plist.items(), key=operator.itemgetter(1, 0), reverse=True)

            # There is a small mistake in the algorithm,
            # it should search not for ALL peaks
            # but only for ALL peaks except itself; position i can not be a noise product of i itself

            n = range(len(psorted))

            for i in n:
                if(psorted[i] is not None):
                    item = psorted[i]
                    for j in n:							# Can be limited to size and -size of self.self.settings.parametersmatrix
                        if((psorted[j] is not None) and (j != i)):
                            item2 = psorted[j]
                            diff = item2[0] - item[0]
                            if diff in self.settings.parameters.matrix:
                                perc = self.settings.parameters.matrix[diff] / 100.0
                                if((perc * item[1]) > item2[1]):
                                    psorted[j] = None

            return {x[0]: x[1] for x in psorted if x is not None}

        def step_04__assemble_fragments(pstart, pstop, pexpectedStart, pexpectedStop):
            """Assemble by peak reconstruction / traceback
            """
            logging.debug("Assembling fragments")

            if len(pstart) >= len(pstop):									# More start than stop positions
                pstopSorted = sort_frequency_dict(pstop)
                for itema in pstopSorted:
                    pos = itema[0]
                    diff = pexpectedStop[pos]
                    predictedPos = pos + diff + 1							# 149 - 50 = 99; 149- 50 + 1 = 100 (example of read aligned to 100,149 (size=50)

                    highest_scoring_position = (0, -1, -1, -1, -1, 99999)

                    for item in sorted(pstart.keys()):
                        if item >= predictedPos - self.settings.parameters.left_padding and item <= predictedPos + self.settings.parameters.right_padding:
                            distance = abs(predictedPos - item)
                            penalty = 1.0 - (distance * 0.09)

                            score = pstart[item] * penalty
                            if (score > highest_scoring_position[0]) or (score == highest_scoring_position[0] and distance < highest_scoring_position[5]):
                                highest_scoring_position = (score, item, pos, pstart[item], pstop[pos], distance)

                    if highest_scoring_position[0] > 0.0:
                        del(pstart[highest_scoring_position[1]])
                        yield ncRNAFragment(highest_scoring_position[1], highest_scoring_position[2], highest_scoring_position[3], highest_scoring_position[4])
            else:															# More stop than start positions
                pstartSorted = sort_frequency_dict(pstart)
                for itema in pstartSorted:
                    pos = itema[0]
                    diff = pexpectedStart[pos]
                    predictedPos = pos + diff  # @todo figure out if this requires << + 1

                    highest_scoring_position = (0, -1, -1, -1, -1, 99999)

                    for item in sorted(pstop.keys(), reverse=True):
                        if item >= predictedPos - self.settings.parameters.left_padding and item <= predictedPos + self.settings.parameters.right_padding:
                            distance = abs(predictedPos - item)
                            penalty = 1.0 - (distance * 0.09)

                            score = pstop[item] * penalty
                            if (score > highest_scoring_position[0]) or (score == highest_scoring_position[0] and distance < highest_scoring_position[5]):
                                highest_scoring_position = (score, pos, item, pstart[pos], pstop[item], distance)  # (Highest score -- a bug... should be 'score', Start, Stop, Reads on start, Reads on stop)

                    if highest_scoring_position[0] > 0.0:
                        del(pstop[highest_scoring_position[2]])
                        yield ncRNAFragment(highest_scoring_position[1], highest_scoring_position[2], highest_scoring_position[3], highest_scoring_position[4])

        # Acquire statistics
        start_positions, stop_positions, start_avg_lengths, stop_avg_lengths = step_01__parse_stats()
        
        print (start_positions)
        print (stop_positions)
        print ("")
        print (start_avg_lengths)
        print (stop_avg_lengths)
        print ("")


        # Finds peaks
        start_positions = step_02__find_peaks(start_positions + [0])
        stop_positions = step_02__find_peaks(stop_positions + [0])

        # Correct / filter noisy peaks
        start_positions = step_03__smooth_filter_peaks(start_positions)
        stop_positions = step_03__smooth_filter_peaks(stop_positions)

        # Trace start and stop positions together and obtain actual peaks
        for fragment in step_04__assemble_fragments(start_positions, stop_positions, start_avg_lengths, stop_avg_lengths):
            yield fragment

    def __iter__(self):
        for fragment in self.predict_fragments():
            yield fragment
