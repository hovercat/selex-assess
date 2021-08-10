#!/usr/bin/env nextflow
"""
========================================================
Groovy Helper Functions
========================================================    
"""
def filter_selex_rounds(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return true;
    // // building regex to check files if they match the rounds specified in YOUR_SELEX
    // round_regex = "(" + params.round_order.join('|') + ")" + params.trim_delimiter + ".*"; // TODO replace trim_delimiter with autodetect
    // return s.matches(round_regex);
    
    return params.round_order.contains(round_name);
}

def get_round_id(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return 0;
    else return params.round_order.indexOf(round_name);
}

"""
========================================================
Make output directory if it doesn't exist
========================================================    
"""
dir_output = file(params.output_dir)
if (!dir_output.exists()) {
    if (!dir_output.mkdir()) {
        println("Couldn't create output directory.")
    }
}


"""
========================================================
Read FASTA-files from specified input directory and order by round_order (if present)
========================================================    
"""
fasta_files = Channel.fromPath(params.input_dir + "/" + params.fasta_pattern, checkIfExists:true, type: "file")
fasta_files
    .map { it -> [get_round_id(it.getSimpleName()), it.getSimpleName(), it].flatten() }
    .filter { filter_selex_rounds(it[1]) }
    .into { fasta_files_filtered; fasta_files_nt  }
   
fasta_files_filtered
    .toSortedList( {a, b -> a[0] <=> b[0] } )
    .collect { it -> return it.collect { it[2] } }
    .set { fasta_files_sorted }

"""
========================================================
Analysing SELEX Enrichment
========================================================
"""

process dereplicate_rpm {
    publishDir "${params.output_dir}/",
        pattern: '*{csv,fasta}',
        mode: "copy"
                
    input:
        file(fasta) from fasta_files_sorted
    output:
        tuple file("aptamers.fasta"), file("aptamers.csv"), file("aptamers.rpm.csv") into selex_dereplicated
	file("aptamers.csv") into selex_top_n
        
    """
        selex_dereplicate_fasta.py -o aptamers.fasta -c aptamers.csv ${fasta}
        selex_rpm.r -i aptamers.csv -o aptamers.rpm.csv
    """
}

process selex_top_n {
    publishDir "${params.output_dir}/",
        pattern: '*.xlsx',
        mode: 'copy'
    input:
        file aptamers_csv from selex_top_n
    output:
        file "top_${params.top_n}.xlsx" into selex_top_n_out
    script:
    """
		selex_top_n.R -i $aptamers_csv -o top_${params.top_n}.xlsx -n ${params.top_n}
    """
}

process assess_selex_enrichment {
    publishDir "${params.output_dir}/analysis.enrichment",
        pattern: '*',
        mode: "copy"
                
    input:
    	tuple file("aptamers.fasta"), file("aptamers.csv"), file("aptamers.rpm.csv") from selex_dereplicated
    output:
        file("enrichment*.csv") into selex_enrichment
        file("enrichment*.png") into selex_enrichment_png
    """        
        # Analyse Log Duplicates
        selex_analyse_log_duplicates.r -i aptamers.csv --log-base 2 --out-unique-csv enrichment.unique.log2.csv --out-csv enrichment.log2.csv
        selex_analyse_log_duplicates.r -i aptamers.csv --log-base 10 --out-unique-csv enrichment.unique.log10.csv --out-csv enrichment.log10.csv
        selex_analyse_log_duplicates.r -i aptamers.rpm.csv --log-base 2 --out-unique-csv enrichment.unique.log2.rpm.csv --out-csv enrichment.log2.rpm.csv
        selex_analyse_log_duplicates.r -i aptamers.rpm.csv --log-base 10 --out-unique-csv enrichment.unique.log10.rpm.csv --out-csv enrichment.log10.rpm.csv
    """
}


"""
========================================================
Analysing Nucleotide Distribution
========================================================
"""
process analyse_selex_nt_distribution {
    publishDir "${params.output_dir}/analysis.nt_distribution",
        pattern: '*.csv',
        mode: "copy"
    input:
        tuple val(round_id), val(round_name), file(fasta) from fasta_files_nt
    output:
        tuple val(round_id), file("${round_name}.nt_distribution.csv") into nt_distribution_round_csv
    script:
    """
        selex_nt_composition.py -i $fasta -o ${round_name}.nt_distribution.csv --DNA -n ${params.random_region}
    """
} 

nt_distr_round_csv_sorted = nt_distribution_round_csv
    .toSortedList( {a -> a[0] } )
    .transpose()
    .last()
    .collect()

process plot_selex_nt_distribution {
    publishDir "${params.output_dir}/analysis.nt_distribution",
        pattern: '*{png,html}',
        mode: "copy"
    input:
        file(round_csv) from nt_distr_round_csv_sorted
    output:
        file("nucleotide_composition.html") into nt_distribution_round_html
        file("*.png") into nt_distribution_round_png
    script:
    """
        selex_nt_composition_plot.r -i $round_csv -o ./nucleotide_composition.html
    """
}
