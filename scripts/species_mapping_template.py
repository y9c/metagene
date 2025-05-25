# Species Mapping Configuration Template
# Add this to your metagene package after generating Parquet files

BUILTIN_SPECIES = {
    # Human
    'GRCh38': 'GRCh38.parquet',
    'hg38': 'hg38.parquet',
    'GRCh37': 'GRCh37.parquet', 
    'hg19': 'hg19.parquet',
    
    # Mouse
    'GRCm39': 'GRCm39.parquet',
    'mm39': 'mm39.parquet',
    'GRCm38': 'GRCm38.parquet',
    'mm10': 'mm10.parquet',
    'mm9': 'mm9.parquet',
    'NCBIM37': 'NCBIM37.parquet',
    
    # Model organisms
    'TAIR10': 'TAIR10.parquet',  # Arabidopsis
    'IRGSP-1.0': 'IRGSP-1.0.parquet',  # Rice
    'WBcel235': 'WBcel235.parquet',  # C. elegans
    'ce11': 'ce11.parquet',  # C. elegans
    'BDGP6.32': 'BDGP6.32.parquet',  # Drosophila
    'dm6': 'dm6.parquet',  # Drosophila
    'GRCz11': 'GRCz11.parquet',  # Zebrafish
    'GRCz10': 'GRCz10.parquet',  # Zebrafish
    'bGalGal1': 'bGalGal1.parquet',  # Chicken
    'Glycine_max_v2.1': 'Glycine_max_v2.1.parquet',  # Soybean
    'R64-1-1': 'R64-1-1.parquet',  # S. cerevisiae
    'sacCer3': 'sacCer3.parquet',  # S. cerevisiae
    'ASM294v2': 'ASM294v2.parquet',  # S. pombe
}
