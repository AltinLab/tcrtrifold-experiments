nextflow.enable.dsl = 2
params.pdb_dir  = null     
params.out_dir  = null    
params.csv_file = null    

if( [params.pdb_dir, params.out_dir, params.csv_file].any{ it == null } )
    exit 1, """
    Please supply:
      --pdb_dir   <dir with PDB files>
      --out_dir   <output dir>
      --csv_file  <CSV metadata>
    """

Channel
    .fromPath( params.csv_file )
    .splitCsv(header:true)
    .map { row ->                            
        def pdb_id     = row.pdbid.trim()
        def organism   = row.organism.trim().toLowerCase()
        def mhc_class  = row.'mhc_class'.toString().trim()  
        tuple( pdb_id, organism, mhc_class )
    }
    .set { meta_ch }

process ParsePDB {
    tag { "${pdb_id}" }
    errorStrategy "ignore"

    input:
    tuple val(pdb_id), val(org), val(mhc_cls)

    script:
    """
    mkdir -p ${params.out_dir}

    . ~/miniconda3/etc/profile.d/conda.sh
    conda activate tcrdock_test

    pdbfile=""
    # a) pattern  <ID>*.pdb   (first match wins)
    for f in ${params.pdb_dir}/${pdb_id}*.pdb ; do
        [[ -f "\$f" ]] && { pdbfile="\$f"; break; }
    done

    python /tgen_labs/altin/alphafold3/runs/tcrtrifold-experiments/data/TCRdock/parse_tcr_pmhc_pdbfile.py \
        --pdbfiles   "\$pdbfile" \
        --organism   ${org} \
        --mhc_class  ${mhc_cls} \
        --out_tsvfile ${params.out_dir}/${pdb_id}.tsv
    """
}

workflow {
    ParsePDB(meta_ch)
}
