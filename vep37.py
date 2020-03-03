import requests, sys
import pprint
import TSO500
import utils
from annotate import add_ckb_variant_info, add_ckb_gene_info, add_clinvar, add_hot_spot_info, is_near_GOF_LOF_mutation, \
    save_annotated_snv, get_pos_from_pdot, get_annotated_snv


def main():
    client = utils.get_mongo_client()
    db = utils.get_database(client,'omni')
    path = 'data/TSO500_UniqueVariants_runs1-6.csv'
    # path = 'data/TSO500_UniqueVariants_4.csv'
    var_list = TSO500.read_tso_unique_variants(path)
    variants_to_analyze = []
    for index, variant in enumerate(var_list):
        key = get_key_from_variant(variant)
        annotated = get_annotated_snv(key,db)
        if annotated == None:
            variants_to_analyze.append(variant)
            if len(variants_to_analyze) >= 100:
                report_on_variants(db, variants_to_analyze)
                variants_to_analyze = []
    if len(variants_to_analyze)>0:
        report_on_variants(db, variants_to_analyze)


def get_best_t_c(t_c,gene,cdot):
    best = t_c[0]
    for conseq in t_c:
        cdot_pos = 'xx'
        if 'cdna_start' in conseq:
            cdot_pos = str(get_pos_from_pdot(cdot_from_consequence(conseq)))

        if conseq['gene_symbol']==gene and cdot==cdot_from_consequence(conseq):
            best = conseq
            break
        elif conseq['gene_symbol']==gene and cdot_pos in cdot:
            best = conseq
            break
    return best


def get_key_from_variant(v):
    key = ''
    pos = int(v['pos_hg19'])
    key += v['chrom']
    key += '_'
    ref = v['ref']
    alt = v['alt']
    if len(ref) == 0:
        # insertion, so: An insertion (of any size) is indicated by start coordinate = end coordinate + 1
        key += str(pos + 1)
        key += '_'
        ref = '-'
    elif len(alt) == 0:
        # deletion: A deletion is indicated by the exact nucleotide coordinates
        l = len(ref) - 1
        key += str(pos)
        key += '_'
        alt = '-'
    elif len(ref) > 1 and len(alt) > 1:
        # deletion: A deletion is indicated by the exact nucleotide coordinates
        key += str(pos)
        key += '_'
    else:
        # variant, sp start and end are equal
        key += str(pos)
        key += '_'

    key += ref
    key += '/'
    key += alt

    return key

def report_on_variants(db, var_list):
    vep_variants = vep_variants_array_from_variants(var_list)
    transcript_consequences = call_vep37(vep_variants)
    for index, variant in enumerate(var_list):
        key = get_key_from_variant(variant)
        if key in transcript_consequences:
            t_c = transcript_consequences[key]
            # best_t_c = t_c[0]
            best_t_c = get_best_t_c(t_c,variant['HGNC_Symbol'],variant['omni_c_dot'])
            parse_one_transcript_consequences(best_t_c, variant, db)

            add_clinvar(variant, db)
            add_hot_spot_info(variant, db)
            is_near_GOF_LOF_mutation(variant, db)
            save_annotated_snv(variant, key, db)
            star = ''
            if variant['HGNC_Symbol'] != variant['gene'] or get_pos_from_pdot(variant['omni_p_dot']) != get_pos_from_pdot(variant['pdot']):
                star = "*"
            print(star,variant['HGNC_Symbol'], variant['gene'], variant['omni_c_dot'], variant['cdot'], variant['omni_p_dot'],
                  variant['pdot'], variant['gene_category'], variant['mutation_type'], variant['in_clinvar'],
                  variant['clinvar_significance'], variant['is_near_GOF_LOF_mutation'])
        else:
            print('missing:',variant['HGNC_Symbol'], variant['omni_c_dot'], variant['omni_p_dot'])


def vep_variants_array_from_variants(var_list):
    variants = []
    for v in var_list:
        variants.append(get_vep_var_as_string(v['chrom'],int(v['pos_hg19']) , v['ref'], v['alt']))
    return variants



def get_vep_var_as_string(chr, pos, ref, alt):
    s = '"'
    s += chr
    s +=  ' '
    if len(ref)==0:
        # insertion, so: An insertion (of any size) is indicated by start coordinate = end coordinate + 1
        s += str(pos+1)
        s += ' '
        s += str(pos)
        s += ' '
        ref = '-'
    elif len(alt)==0:
        # deletion: A deletion is indicated by the exact nucleotide coordinates
        l = len(ref) - 1
        s += str(pos)
        s += ' '
        s += str(pos+l)
        s += ' '
        alt= '-'
    elif len(ref) > 1 and len(alt) > 1:
        # deletion: A deletion is indicated by the exact nucleotide coordinates
        l = len(ref) - 1
        s += str(pos)
        s += ' '
        s += str(pos + l)
        s += ' '
    else:
        # variant, sp start and end are equal
        s += str(pos)
        s += ' '
        s += str(pos)
        s += ' '

    s += ref
    s += '/'
    s += alt
    s += '"'
    return s


def get_pdot_from_hgvsp(hgsvp):
    pdot = hgsvp[hgsvp.find(':p')+3:]
    pdot = pdot.replace('Ala','A')
    pdot = pdot.replace('Cys','C')
    pdot = pdot.replace('Asp','D')
    pdot = pdot.replace('Glu','E')
    pdot = pdot.replace('Phe','F')
    pdot = pdot.replace('Gly','G')
    pdot = pdot.replace('His','H')
    pdot = pdot.replace('Ile','I')
    pdot = pdot.replace('Lys','K')
    pdot = pdot.replace('Leu','L')
    pdot = pdot.replace('Met','M')
    pdot = pdot.replace('Asn','N')
    pdot = pdot.replace('Pro','P')
    pdot = pdot.replace('Gln','Q')
    pdot = pdot.replace('Arg','R')
    pdot = pdot.replace('Ser','S')
    pdot = pdot.replace('Thr','T')
    pdot = pdot.replace('Val','V')
    pdot = pdot.replace('Trp','W')
    pdot = pdot.replace('Tyr','Y')
    pdot = pdot.replace('Ter','*')
    return pdot


def parse_one_transcript_consequences(consequence,variant,db):
    protein_altering_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                                 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
                                 'transcript_amplification',
                                 'inframe_insertion', 'inframe_deletion', 'missense_variant',
                                 'protein_altering_variant']
    truncating_variants = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
                           'frameshift_variant', 'start_lost']

    variant['gene'] = consequence['gene_symbol']
    add_ckb_gene_info(variant, db)
    # add_uniprot_info(variant)
    mutation_type = consequence['consequence_terms'][0]
    cdot = cdot_from_consequence(consequence)
    variant['cdot'] = cdot

    variant['pdot'] = ''
    if 'protein_start' in consequence:
        protein_start = consequence['protein_start']
        variant['protein_start'] = protein_start
        # add_region_hit(variant)
        # add_hot_spot_info(variant)
        if 'hgvsp' in consequence:
            variant['pdot'] = get_pdot_from_hgvsp(consequence['hgvsp'])
        # add_ckb_variant_info(variant, db)
        variant['full_name'] = variant['gene'] + ' ' + variant['pdot']

    variant['mutation_type'] = mutation_type
    variant['is_protein_altering'] = mutation_type in protein_altering_variants
    variant['is_truncating_variants'] = mutation_type in truncating_variants
    if 'polyphen_prediction' in consequence:
        variant['polyphen_prediction'] = consequence['polyphen_prediction']
    if 'sift_prediction' in consequence:
        variant['sift_prediction'] = consequence['sift_prediction']
    variant['predicted_deleterious'] = False
    if 'polyphen_prediction' in variant and 'sift_prediction' in variant:
        if 'damaging' in variant['polyphen_prediction'] and 'deleterious' in variant[
            'sift_prediction']:
            variant['predicted_deleterious'] = True

    add_ckb_variant_info(variant, db)


def cdot_from_consequence(consequence):
    cdot = '-'
    if 'hgvsc' in consequence:
        cdot = consequence['hgvsc'][consequence['hgvsc'].find(':') + 1:]
    return cdot


def call_vep37(variants):
    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    body = '{ "variants" : ['
    v_body = ''
    for v in variants:
        if len(v_body)>0:
            v_body += ', '
        v_body += v
    body += v_body + '] , "canonical": "True", "hgvs":"True" }'
    # print(body)
    r = requests.post(server + ext, headers=headers,data=body)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    transcript_consequences = {}
    for elem in decoded:
        if 'transcript_consequences' in elem:
            key = elem['id']
            transcript_consequences[key] = elem['transcript_consequences']

    return transcript_consequences


if __name__ == "__main__":
    main()