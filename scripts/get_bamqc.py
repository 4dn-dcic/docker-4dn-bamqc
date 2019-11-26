import pysam
import json
import click
import collections
import utilities


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
    pairs_tags_dict = collections.OrderedDict()
    pairs_tags_dict['Minor Contigs'] = 0
    pairs_tag_list = []
    main_chroms = [x for x in range(len(chrmsizes))] # the chromsizes file
    candidates = []

    for read in bamfile:
        candidate_info = {'id': read.query_name,
                          'chrm': read.reference_id,
                          'pos': read.pos,
                          'chrm_next': read.next_reference_id,
                          'pos_next': read.pnext,
                          'tag': read.get_tag('Yt'),
                          'chrm_pos': (read.reference_id, read.pos),
                          'chrm_pos_next': (read.next_reference_id, read.pnext)}

        if candidate_info['tag'] not in pairs_tags_dict.keys():
            pairs_tags_dict[candidate_info['tag']] = 0

        if prev_read_id is not None and read.query_name != prev_read_id:
            # This is a group of lines with the same id. Send to process
            utilities.process_read_ids_block(candidates, pairs_tag_list, pairs_tags_dict, main_chroms)
            pairs_tag_list = []
            candidates = []

        pairs_tag_list.append(candidate_info['tag'])
        prev_read_id = read.query_name
        candidates.append(candidate_info)

    if len(pairs_tag_list) != 0:
        utilities.process_read_ids_block(candidates, pairs_tag_list, pairs_tags_dict, main_chroms)

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
