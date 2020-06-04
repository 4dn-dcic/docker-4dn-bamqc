import pysam
import json
import click
import collections


@click.command()
@click.argument('infile')
@click.argument('chromsizes')
@click.argument('outdir')
@click.argument('filename')
def main(infile, chromsizes, outdir, filename):
    f = infile
    outpath = outdir + '/' + filename + '.json'
    bamfile = pysam.AlignmentFile(f, "rb")
    with open(chromsizes, 'r') as opf:
        chrmsizes = opf.readlines()

    prev_read_id = None
    pairs_tags_count_dict = collections.OrderedDict()
    pairs_tags_count_dict['Minor Contigs'] = 0 # Keeps track of minor contigs for UU, RU and UR flag types
    pairs_tag_list = []
    main_chroms = [a_line.split('\t')[0] for a_line in chrmsizes]  # the chromsizes file
    lines_same_id = []

    for bamline in bamfile:
        # Getting bamline information
        bamline_info = {'id': bamline.query_name,
                        'chrm': bamline.reference_name,
                        'tag': bamline.get_tag('Yt'),
                        'secondary_alignment': bamline.is_secondary}

        if bamline_info['tag'] not in pairs_tags_count_dict.keys():
            pairs_tags_count_dict[bamline_info['tag']] = 0

        if prev_read_id is not None and bamline.query_name != prev_read_id:
            # A line with a new read_id is being reached
            # Group lines with the previous read_id into reads and count number of reads per tag type.
            lines_per_tag_groups = []  # grouping lines per tag type
            for tag in set(pairs_tag_list):
                lines_per_tag_group = []
                for line_same_id in lines_same_id:
                    if line_same_id['tag'] == tag:
                        lines_per_tag_group.append(line_same_id)
                lines_per_tag_groups.append(lines_per_tag_group)

            for tag_group in lines_per_tag_groups:  # Processing each group of tag types
                read_num = 0
                track_minor_contigs = False
                subset_value = 0 # how many lines correspond to a single read for each type of flag, use this info to keep track of minor contigs
                if tag_group[0]['tag'] == "UU":  # count minor contigs for this flag type
                    track_minor_contigs = True
                    subset_counter = 2
                    subset_value = 2  # a read that has the UU flag is represented by 2 lines only.
                    is_minor_contig = False
                elif (tag_group[0]['tag'] == "RU" or tag_group[0]['tag'] == "UR"):  # count minor contigs for this flag type
                    track_minor_contigs = True
                    subset_value = 3  # a read that has the RU and UR flag is represented by 3 lines only.
                    subset_counter = 3
                    is_minor_contig = False

                for i, a_line in enumerate(tag_group):
                    if not a_line['secondary_alignment']: # a read has only 2 primary alignments, the rest are secondary.
                        read_num = read_num + 1

                    # for UU, UR and RU flags check if a line within a read is mapped to a non main chromosome,
                    # this count will be added to minor contigs count and removed from filtered reads count in the final output
                    if track_minor_contigs:

                        if a_line['chrm'] not in main_chroms:
                            is_minor_contig = True

                        if (i + 1) == subset_counter:
                            if is_minor_contig:
                                pairs_tags_count_dict['Minor Contigs'] = pairs_tags_count_dict['Minor Contigs'] + 1
                                is_minor_contig = False
                            subset_counter = subset_counter + subset_value

                num_reads_tag_group = read_num // 2 # divide the number of primary alignments by 2 to obtain the number of reads per tag type
                pairs_tags_count_dict[tag_group[0]['tag']] = pairs_tags_count_dict[tag_group[0]['tag']] + num_reads_tag_group

            pairs_tag_list = []
            lines_same_id = []

        pairs_tag_list.append(bamline_info['tag'])
        prev_read_id = bamline.query_name
        lines_same_id.append(bamline_info)

    if len(pairs_tag_list) != 0:  # There are some lines left to process, repeat processing as above
        lines_per_tag_groups = []
        for tag in set(pairs_tag_list):
            lines_per_tag_group = []
            for line_same_id in lines_same_id:
                if line_same_id['tag'] == tag:
                    lines_per_tag_group.append(line_same_id)
            lines_per_tag_groups.append(lines_per_tag_group)

        for tag_group in lines_per_tag_groups:
            read_num = 0
            track_minor_contigs = False
            subset_value = 0
            if tag_group[0]['tag'] == "UU":
                track_minor_contigs = True
                subset_counter = 2
                subset_value = 2
                is_minor_contig = False
            elif (tag_group[0]['tag'] == "RU" or tag_group[0]['tag'] == "UR"):
                track_minor_contigs = True
                subset_value = 3
                subset_counter = 3
                is_minor_contig = False

            for i, a_line in enumerate(tag_group):
                if not a_line['secondary_alignment']:
                    read_num = read_num + 1
                if track_minor_contigs:

                    if a_line['chrm'] not in main_chroms:
                        is_minor_contig = True

                    if (i + 1) == subset_counter:
                        if is_minor_contig:
                            pairs_tags_count_dict['Minor Contigs'] = pairs_tags_count_dict['Minor Contigs'] + 1
                            is_minor_contig = False
                        subset_counter = subset_counter + subset_value

            num_reads_tag_group = read_num // 2
            pairs_tags_count_dict[tag_group[0]['tag']] = pairs_tags_count_dict[tag_group[0]['tag']] + num_reads_tag_group

    # QC Report
    report_dict = collections.OrderedDict()
    total_reads = 0
    rescued_chimeras = 0
    unmmapped = 0
    multi = 0
    duplicates = 0
    unique = 0

    for key, val in pairs_tags_count_dict.items(): # Categorize each flag type
        if key != "Minor Contigs":
            total_reads = total_reads + val
        if key == "RU" or key == "UR":
            rescued_chimeras = rescued_chimeras + val
        elif "N" in key or "X" in key:
            unmmapped = unmmapped + val
        elif key == "DD":
            duplicates = val
        elif "M" in key and key != "NM":
            multi = multi + val
        elif key == "UU":
            unique = val

    filtered = unique + rescued_chimeras - pairs_tags_count_dict['Minor Contigs']

    report_dict['Total Reads'] = total_reads
    report_dict['% of reads mapped'] = (1 - (unmmapped / total_reads)) * 100
    report_dict['% of duplicates'] = duplicates/total_reads * 100
    report_dict['% of reads filtered'] = filtered/total_reads * 100
    report_dict['% of rescued chimeras'] = rescued_chimeras/total_reads * 100
    report_dict['Filtered Reads'] = filtered
    report_dict.update(pairs_tags_count_dict)

    with open(outpath, 'w') as outfile:
        json.dump(report_dict, outfile, indent=2)


if __name__ == "__main__":
    main()
