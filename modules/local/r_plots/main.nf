process R_PLOT {
    tag "tree_amr_plot"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path tree
    path metadata
    path rgi_txts

    output:
    path "plot_tree.pdf", emit: tree_pdf
    path "plot_tree.png", emit: tree_png
    path "plot_AMR.pdf",  emit: amr_pdf
    path "plot_AMR.png",  emit: amr_png

    script:
    """
    Rscript ${moduleDir}/analysis.R ${tree} ${metadata} plot ${rgi_txts}
    """
}