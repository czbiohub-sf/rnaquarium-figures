def parse_gtf_attributes(attribute_string):
    """
    Parse the GTF attribute string and return a dictionary of key-value pairs.
    Example of attribute_string:
       'gene_id "ENSG00000198947"; transcript_id "ENST00000382353"; ...'
    """
    attr_dict = {}
    # Split on semicolon to get each key-value pair
    for attribute in attribute_string.strip().split(';'):
        attribute = attribute.strip()
        if not attribute:
            # Skip empty parts
            continue

        # Each attribute typically looks like key "value"
        # We can split on the first space to separate the key from the quoted value
        parts = attribute.split(' ', 1)
        if len(parts) != 2:
            continue

        key, value = parts
        # Remove surrounding quotes from the value
        value = value.strip('"')
        attr_dict[key] = value

    return attr_dict


def extract_gene_biotypes(gtf_file):
    """
    Read a GTF file and return a dictionary mapping gene_id -> gene_biotype.
    """
    gene_biotypes = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            # Split the GTF line into the 9 columns
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            feature_type = columns[2]  # e.g. "gene", "transcript", "exon", etc.
            attributes_str = columns[8]

            # We only want to look at lines describing a gene
            if feature_type != 'gene':
                continue

            # Parse the attributes into a dictionary
            attr_dict = parse_gtf_attributes(attributes_str)

            # Extract gene_id and gene_biotype if present
            if 'gene_id' in attr_dict:
                t_id = attr_dict['gene_id']
                # Some lines may not have gene_biotype, so use get() with a default
                biotype = attr_dict.get('gene_biotype', 'N/A')
                gene_biotypes[t_id] = biotype

    return gene_biotypes