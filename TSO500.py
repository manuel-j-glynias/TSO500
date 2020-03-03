import annotate
import utils
import csv
import numbers



def get_variant_from_row(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['omni_c_dot'] = row['c_dot']
    variant['omni_p_dot'] = row['p_dot']
    variant['omni_MutationType'] = row['MutationType']
    chrom = row['Chr']
    if isinstance(chrom, numbers.Number):
        chrom = int(chrom)

    variant['chrom'] = str(chrom)
    variant['ref'] = row['Ref']
    variant['alt'] = row['Alt']
    # old pos
    variant['pos_hg19'] = row['Position']
    if not isinstance(variant['pos_hg19'], numbers.Number):
        variant['pos_hg19'] = variant['pos_hg19'].replace(u'\ufeff', '')
    variant['type'] = 'snv'
    variant['report_status'] = 'not_reported'
    return variant


def read_tso_unique_variants(path):
    variants = []
    with open(path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            variants.append(get_variant_from_row(row))
    return variants


