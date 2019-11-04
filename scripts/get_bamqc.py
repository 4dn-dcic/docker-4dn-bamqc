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

    prev_read = None
    pairs_tags_dict = collections.OrderedDict()
    pairs_tags_dict['Minor Contigs'] = 0
    pairs_tag_list = []
    main_chroms = [x for x in range(len(chrmsizes))]
    prev_ignored_read = False

    for read in bamfile:
        pairs_tag = read.get_tag('Yt')
        chrm = read.reference_id
        if pairs_tag not in pairs_tags_dict.keys():
            pairs_tags_dict[pairs_tag] = 0

        if prev_read is not None and read.query_name != prev_read:
            # This is a new read id, Add the information of the previous read id
            assert (len(set(pairs_tag_list)) == 1), "Read %s has more than one pairs flag " % prev_read
            pairs_tags_dict[pairs_tag_list[0]] = pairs_tags_dict[pairs_tag_list[0]] + 1

            if prev_ignored_read:
                pairs_tags_dict['Minor Contigs'] = pairs_tags_dict['Minor Contigs'] + 1
            pairs_tag_list = []
            prev_ignored_read = False

        pairs_tag_list.append(pairs_tag)
        prev_read = read.query_name
        if (pairs_tag == "UU" or pairs_tag == "RU" or pairs_tag == "UR") and (chrm not in main_chroms):
            prev_ignored_read = True

    if len(pairs_tag_list) != 0:
        assert (len(set(pairs_tag_list)) == 1), "Read %s has more than one pairs flag " % prev_read
        pairs_tags_dict[pairs_tag_list[0]] = pairs_tags_dict[pairs_tag_list[0]] + 1
        if prev_ignored_read:
            pairs_tags_dict['Minor Contigs'] = pairs_tags_dict['Minor Contigs'] + 1

    # QC Report
    report_dict = collections.OrderedDict()
    total_reads = 0
    rescued_chimeras = 0
    unmmapped = 0
    multi = 0
    duplicates = 0
    unique = 0

    for key, val in pairs_tags_dict.items():
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

    filtered = unique + rescued_chimeras - pairs_tags_dict['Minor Contigs']

    report_dict['Total Reads'] = total_reads
    report_dict['% of reads mapped'] = (1 - (unmmapped / total_reads)) * 100
    report_dict['% of duplicates'] = duplicates/total_reads * 100
    report_dict['% of reads filtered'] = filtered/total_reads * 100
    report_dict['% of rescued chimeras'] = rescued_chimeras/total_reads * 100
    report_dict['Filtered Reads'] = filtered
    report_dict.update(pairs_tags_dict)

    with open(outpath, 'w') as outfile:
        json.dump(report_dict, outfile, indent=2)


if __name__ == "__main__":
    main()
