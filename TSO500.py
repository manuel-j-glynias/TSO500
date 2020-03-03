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


def annotate_TSO_snv(strands,variant,db):
    annotated_snv = annotate.get_annotated_snv(variant,db)
    if annotated_snv is None:
        annotate.handle_liftover(strands, variant)
        annotate.call_annotate(variant, db)
        if 'pdot' not in variant:
            variant['pdot'] = ''
        annotate.save_annotated_snv(variant,db)
        annotated_snv = variant
    return annotated_snv


def read_strands(path):
    strands = {}
    with open(path) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            strands[row[0]] = row[1]
    return strands

def main():
    client = utils.get_mongo_client()
    db = utils.get_database(client,'omni')

    db.drop_collection('ann_var')
    db.create_collection('ann_var')

    strands = read_strands('data/strands.tsv')
    variants =read_tso_unique_variants()
    for index, variant in enumerate(variants):
        print(variant['chrom'],variant['pos_hg19'],variant['ref'],variant['alt'])
        annotated_variant = annotate_TSO_snv(strands, variant, db)
        print(index,annotated_variant)
        print(annotated_variant['pdot'], annotated_variant['gene_category'],annotated_variant['in_clinvar'],annotated_variant['clinvar_significance'],annotated_variant['is_near_GOF_LOF_mutation'])
        print()

if __name__ == "__main__":
    main()