def map_reads_to_groups(candidates):
    '''
    when there are more than 4 lines of in which they have the same read_id and the same tag.
    Verify if all the lines correspond to the same read or multiple reads
    args: candidates is the group of lines with the same read_id and tag
    returns a list of list representing the groups the lines belong to, each group correspond to a single read.
    '''
    all_groups = []

    # check if the reads have the same position
    candidates_map = {}

    for idx, candidate in enumerate(candidates):
        assert (candidate['chrm_pos'] not in candidates_map), "There are multiple lines with same read id (%s), tag (%s) and position (%s)" % (str(candidate['id']), str(candidate['tag']), str(candidate['chrm_pos']))
        candidates_map[candidate['chrm_pos']] = idx

    real_group_indx = []
    candidate_group_indx = []
    belong = False
    for candidate in candidates:
        current_idx = candidates_map[candidate['chrm_pos']]
        if current_idx not in real_group_indx:
            if not real_group_indx:
                real_group_indx.append(current_idx)
            belong = False
            next_idx = candidates_map[candidate['chrm_pos_next']]
            candidate_group_indx.append(current_idx)
            while next_idx not in candidate_group_indx:
                candidate_group_indx.append(next_idx)
                next_idx = candidates_map[candidates[next_idx]['chrm_pos_next']]

            for index in real_group_indx:
                if index in candidate_group_indx:
                    belong = True
                    real_group_indx = list(set(real_group_indx + candidate_group_indx))
                    candidate_group_indx = []
                    break

        if not belong:
            real_group = [candidates[i] for i in real_group_indx]
            all_groups.append(real_group)

            real_group_indx = [i for i in candidate_group_indx]
            candidate_group_indx = []
    return all_groups


def process_read_ids_block(candidates, pairs_tag_list, pairs_tags_dict, main_chroms):
    '''
        This function process lines with the same read_id, which represents a single read
        and records the pairs flag
        It also handle cases in which multiple lines have the same read_id but correspond to
        different reads.
        args:
            candidates: a list of dicts, each dict is a single line with useful information.
            pairs_tag_list: a list of the tags in candidates
            pairs_tags_dict: dictionary that hold the count for each pairs flag
            main_chroms: a list containing the numeric indices for the main chromosomes
    '''
    tag_groups = []
    for tag in set(pairs_tag_list):
        tag_group = []
        for candidate in candidates:
            if candidate['tag'] == tag:
                tag_group.append(candidate)
        tag_groups.append(tag_group)

    for tag_group in tag_groups:
        read_groups = []
        if len(tag_group) <= 3:  # This is a single read for sure
            read_groups.append(tag_group)
        elif len(tag_group) == 4 and tag_group[0]['tag'] == "WW":
            read_groups.append(tag_group)
        elif (len(tag_group) % 2 == 0 or len(tag_group) % 3 == 0) and (tag_group[0]['tag'] != "DD" and tag_group[0]['tag'] != "WW"):
            if "R" in tag_group[0]['tag']:
                for x in range(0, len(tag_group), 3):
                    read_group = []
                    for i, tag_member in enumerate(tag_group):
                        if i > x and i < x + 3:
                            read_group.append(tag_member)
                    read_groups.append(read_group)
            if "UU" in tag_group[0]['tag']:
                for x in range(0, len(tag_group), 2):
                    read_group = []
                    for i, tag_member in enumerate(tag_group):
                        if i > x and i < x + 2:
                            read_group.append(tag_member)
                    read_groups.append(read_group)
        else:
            read_groups = map_reads_to_groups(tag_group)

        for group in read_groups:
            pairs_tags_dict[group[0]['tag']] = pairs_tags_dict[group[0]['tag']] + 1
            if (group[0]['tag'] == "UU" or group[0]['tag'] == "RU" or group[0]['tag'] == "UR"):
                for member in group:
                    if member['chrm'] not in main_chroms:
                        pairs_tags_dict['Minor Contigs'] = pairs_tags_dict['Minor Contigs'] + 1
                        break
