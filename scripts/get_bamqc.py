import pysam
import json
import click


@click.command()
@click.argument('infile')
@click.argument('outdir')
@click.argument('filename')
def main(infile, outdir, filename):
    f = infile
    outpath = outdir + '/' + filename + '.json'
    bamfile = pysam.AlignmentFile(f, "rb")

    prev_read = None
    pairs_tags_dict = {}  # Or just have the dictionary with all flags (depends on qc)
    pairs_tag_list = []

    for read in bamfile:
        pairs_tag = read.get_tag('Yt')
        if pairs_tag not in pairs_tags_dict.keys():
            pairs_tags_dict[pairs_tag] = 0

        if prev_read is not None and read.query_name != prev_read:
            # This is a new read id, Add the information of the previous read id
            assert (len(set(pairs_tag_list)) == 1), "Read %s has more than one pairs flag " % prev_read
            pairs_tags_dict[pairs_tag_list[0]] = pairs_tags_dict[pairs_tag_list[0]] + 1
            pairs_tag_list = []

        pairs_tag_list.append(pairs_tag)
        prev_read = read.query_name

    with open(outpath, 'w') as outfile:
        json.dump(pairs_tags_dict, outfile, indent=2)


if __name__ == "__main__":
    main()
