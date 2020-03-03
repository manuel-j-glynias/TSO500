import csv

import requests, sys


def fetch_gene_id_by_gene_name(gene_name):
    requestURL = "http://mygene.info/v3/query?species=human&q=symbol:" + gene_name
    r = requests.get(requestURL, headers={"Accept": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.json()
    # pprint.pprint(responseBody)
    if 'hits' in responseBody and len(responseBody['hits'])>0:
        body = responseBody['hits'][0]['entrezgene']
    else:
        body = None
    return body

def fetch_gene_info_by_gene_id(gene_id):
  requestURL = "http://mygene.info/v3/gene/" + gene_id
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()

  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody


def read_tso_genes():
    genes = []
    with open('data/TSO500_UniqueVariants_runs1-6.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene = row['HGNC_Symbol']
            if gene not in genes:
                genes.append(gene)
    return genes


def create_hgnc_gene_name_dict():
    genes = {}
    genes['C11orf30'] = 'VSIR'
    genes['FAM175A'] = 'ABRAXAS1'
    genes['FAM46C'] = 'TENT5C'
    genes['GPR124'] = 'ADGRA2'
    genes['H3F3C'] = 'H3-5'
    genes['HIST1H1C'] = 'H1-2'
    genes['HIST1H2BB'] = 'H2BC3'
    genes['HIST1H2BD'] = 'H2BC5'
    genes['HIST1H3A'] = 'H3C1'
    genes['HIST1H3E'] = 'H3C6'
    genes['HIST1H3F'] = 'H3C7'
    genes['HIST1H3G'] = 'H3C8'
    genes['HIST2H3D'] = 'H3C13'
    genes['HIST3H3'] = 'H3-4'
    genes['LOC101926927'] = 'withdrawn'
    genes['MARCH9'] = 'MARCHF9'
    genes['MRE11A'] = 'MRE11'
    genes['PAK7'] = 'PAK5'
    genes['PARK2'] = 'PRKN'
    genes['RFWD2'] = 'COP1'
    genes['WISP3'] = 'CCN6'
    return genes

def main():
    gene_dict = create_hgnc_gene_name_dict()
    genes = read_tso_genes()
    for gene in genes:
        gene_name = gene
        if gene in gene_dict:
            gene_name = gene_dict[gene]
        if gene_name=='withdrawn':
            print(gene + '\t' + 'Forward')
        else:
            gene_id = fetch_gene_id_by_gene_name(gene_name)
            if gene_id == None:
                print('*******' +gene_name)
            else:
                gene_info = fetch_gene_info_by_gene_id(gene_id)
                if 'genomic_pos_hg19' in gene_info:
                    genomic_pos = gene_info['genomic_pos_hg19']
                    if isinstance(genomic_pos, list):
                        genomic_pos = genomic_pos[0]
                    strand = "Forward"
                    if genomic_pos['strand']==-1:
                        strand = "Reverse"
                    if gene != gene_name:
                        print(gene_name + '\t' + strand)
                    print(gene + '\t' + strand)
                else:
                    print('#######' + gene)
                    if 'exons_hg19' in gene_info:
                        strand = "Forward"
                        if gene_info['exons_hg19'][0]['strand']==-1:
                            strand = "Reverse"
                        print(gene + '\t' + strand)

if __name__ == "__main__":
    main()