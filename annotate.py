import json
import requests
import re



def add_ckb_gene_info(variant, db):
    gene = variant['gene']
    mycol = db["gene_info"]
    myquery = {'gene': gene}
    mydoc = mycol.find_one(myquery)
    variant['gene_description'] = ''
    variant['gene_category'] = ''
    if mydoc is not None:
        if 'description' in mydoc:
            variant['gene_description'] = mydoc['description']
        if 'category' in mydoc:
            variant['gene_category'] = mydoc['category']

def get_pos_from_pdot(pDot):
    pos = 1
    p = re.compile('\d+\D')
    m = p.search(pDot)
    if m:
        pos = int(m.group()[:-1])
    return pos

def add_ckb_variant_info(variant, db):
    variant['ckb_id'] = ''
    variant['variant_description'] = ''
    variant['protein_effect'] = ''
    variant['variant_type'] = ''
    variant['gDot'] = ''
    variant['cDot'] = ''
    variant['cdot_pos'] = ''
    variant['pdot_pos'] = ''
    variant['gene_id'] = ''
    variant['is_gain_of_function'] = False
    variant['is_loss_of_function'] = False
    if 'pdot' in variant:
        gene = variant['gene']
        pdot = variant['pdot']
        variant['pdot_pos'] = get_pos_from_pdot(pdot)
        full_name = gene + ' ' + pdot
        mycol = db["variant_info"]
        myquery = {'full_name': full_name}
        mydoc = mycol.find_one(myquery)

        if mydoc is not None:
            variant_info_from_ckb(mydoc, variant)
        else:
            if variant['gene_category']=='Tumor Suppressor Gene' and (variant['is_truncating_variants'] or variant[
                'predicted_deleterious']):
                myquery = {'full_name': gene + ' inact mut'}
                mydoc = mycol.find_one(myquery)
                if mydoc is not None:
                    variant_info_from_ckb(mydoc, variant)
            if mydoc is None:
                myquery = {'full_name': gene + ' mutant'}
                mydoc = mycol.find_one(myquery)
                if mydoc is not None:
                    variant_info_from_ckb(mydoc, variant)

def variant_info_from_ckb(mydoc, variant):
    if 'ckb_id' in mydoc:
        variant['ckb_id'] = mydoc['ckb_id']
    if 'description' in mydoc:
        variant['variant_description'] = mydoc['description']
    if 'protein_effect' in mydoc:
        variant['protein_effect'] = mydoc['protein_effect']
        variant['is_gain_of_function'] = 'gain of function' in variant['protein_effect']
        variant['is_loss_of_function'] = 'loss of function' in variant['protein_effect']
    if 'variant_type' in mydoc:
        variant['variant_type'] = mydoc['variant_type']
    if 'gDot' in mydoc and mydoc['gDot'] is not None:
        variant['gDot'] = mydoc['gDot']
    if 'cDot' in mydoc and mydoc['cDot'] is not None:
        variant['cDot'] = mydoc['cDot']
    if 'cdot_pos' in mydoc and mydoc['cdot_pos'] is not None:
        variant['cdot_pos'] = mydoc['cdot_pos']
    if 'pdot_pos' in mydoc and mydoc['pdot_pos'] is not None:
        variant['pdot_pos'] = mydoc['pdot_pos']
    if 'gene_id' in mydoc:
        variant['gene_id'] = mydoc['gene_id']




def get_num_cv_nonbenign_between(gene, begin, end, db):
    mycol = db["clinvar"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               'is_majority_vote_not_benign': True}
    count = mycol.find(myquery).count()
    return count


def get_num_cv_pathogenic_between(gene, begin, end, db):
    mycol = db["clinvar"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               'is_majority_vote_pathogenic': True}
    count = mycol.find(myquery).count()
    return count


def has_cv_pathogenic_or_GOF_LOF_variants_on_both_sides(gene, pos, window, db):
    b = False
    n = get_num_cv_pathogenic_between(gene, pos - window, pos, db)
    if n == 0:
        n = get_num_ckb_gof_lof_variants_between(gene, pos - window, pos, db)
    m = get_num_cv_pathogenic_between(gene, pos, pos + window, db)
    if m == 0:
        m = get_num_ckb_gof_lof_variants_between(gene, pos, pos + window, db)
    if n > 0 and m > 0:
        b = True
    return b, n, m



def get_num_ckb_gof_lof_variants_between(gene, begin, end, db):
    mycol = db["variant_info"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               '$or': [{'protein_effect': 'loss of function'}, {'protein_effect': 'loss of function - predicted'},{'protein_effect': 'gain of function'}, {'protein_effect': 'gain of function - predicted'}]}
    count = mycol.find(myquery).count()
    return count

def get_num_ckb_lof_variants_between(gene, begin, end, db):
    mycol = db["variant_info"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               '$or': [{'protein_effect': 'loss of function'}, {'protein_effect': 'loss of function - predicted'}]}
    count = mycol.find(myquery).count()
    return count


def get_num_ckb_gof_variants_between(gene, begin, end, db):
    mycol = db["variant_info"]
    myquery = {'gene': gene, 'pdot_pos': {'$gte': begin, '$lte': end},
               '$or': [{'protein_effect': 'gain of function'}, {'protein_effect': 'gain of function - predicted'}]}
    count = mycol.find(myquery).count()
    return count


def add_hot_spot_info(variant,db):
    gene = variant['gene']
    mycol = db["hotspots"]
    myquery = {'gene': gene}
    hotspots = []
    if 'pdot_pos' in variant:
        mydocs = mycol.find(myquery).sort("start")
        if mydocs is not None:
            for doc in mydocs:
                if str.isdigit(doc['begin']) and str.isdigit(doc['end']):
                    region = {}
                    region['begin'] = int(doc['begin'])
                    region['end'] = int(doc['end'])
                    if region['begin'] <= variant['pdot_pos'] and region['end'] >= variant['pdot_pos']:
                        hotspots.append(region)
    variant['hotspots'] = hotspots


def add_clinvar(variant, db):
    gene = variant['gene']
    variant['in_clinvar'] = False
    variant['clinvar_significance'] = '-'
    variant['is_clinvar_not_benign'] = '-'
    variant['is_clinvar_benign'] = '-'
    variant['is_clinvar_pathogenic'] = '-'

    if 'pdot' in variant:
        pdot = variant['pdot']
        mycol = db["clinvar"]
        myquery = {'gene': gene, 'pDot': pdot}
        mydoc = mycol.find_one(myquery)
        if (mydoc != None):
            variant['in_clinvar'] = True
            variant['clinvar_significance'] = mydoc['significance']
            variant['clinvar_explain'] = mydoc['explain']
            variant['variant_id'] = mydoc['variant_id']
            variant['is_clinvar_not_benign'] = mydoc['is_majority_vote_not_benign']
            variant['is_clinvar_benign'] = mydoc['is_majority_vote_benign']
            variant['is_clinvar_pathogenic'] = mydoc['is_majority_vote_pathogenic']




def is_near_GOF_LOF_mutation(variant, db):
    window = 10
    gene = variant['gene']
    variant['is_near_GOF_LOF_mutation'] = False
    if 'pdot_pos' in variant:
        pos = variant['pdot_pos']
        if (gene != None and pos != None and pos != ''):
            b, upstream, downstream = has_cv_pathogenic_or_GOF_LOF_variants_on_both_sides(gene, int(pos), window, db)
            variant['is_near_GOF_LOF_mutation'] = b
            variant['num_upstream_GOF_LOF_mutations'] = upstream
            variant['num_downstream_GOF_LOF_mutations'] = downstream






def get_annotated_snv(key, db):
    ann_var = None
    mycol = db["ann_var"]
    # key = get_key_for_variant(variant)
    myquery = {'key': key}
    mydoc = mycol.find_one(myquery)
    if mydoc is not None:
        json_string = mydoc['value']
        ann_var = json.loads(json_string)
    return ann_var


def get_key_for_variant(variant):
    key = variant['chrom'] + '_' + str(variant['pos_hg19']) + '_' + variant['alt']
    return key

def save_annotated_snv(variant, key, db):
#    key = get_key_for_variant(variant)
    json_string = json.dumps(variant)
    mycol = db["ann_var"]
    kv_dict = {'key': key, 'value':json_string}
    mycol.insert_one(kv_dict)

