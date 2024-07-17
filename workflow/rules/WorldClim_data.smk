
rule get_world_clim:
    input: 
    output: temp("../resources/climate_data_{decade}_{variable}_res2.5.zip")
    params: source = "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.06/2.5m/wc2.1_cruts4.06_2.5m_{variable}_{decade}.zip"
    wildcard_constraints: decade = "(1990-1999|2000-2009|2010-2019)", variable = "(tmin|tmax|prec)"
    shell: "wget -O {output} {params.source}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule extract_years:
    input: expand("../resources/climate_data_{decade}_{variable}_res2.5.zip", variable = ["tmin", "tmax", "prec"], decade = ["1990-1999","2000-2009","2010-2019"])
    output: temp(expand("../resources/climate_data/wc2.1_2.5m_{variable}_{year}-{month}.tif", year = ["1997","2009","2010","2017"], month = ["{:02}".format(i) for i in range(1, 13)], variable = ["tmin", "tmax", "prec"]))
    wildcard_constraints: decade = "(1990-1999|2000-2009|2010-2019)", variable = "(tmin|tmax|prec)", year = "(1997|2009|2010|2017)", month = range(1,13)
    shell: "bash scripts/bash/unzip_worldclim.sh {input}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#para rodar essa regra precisa usar com '--use-conda' (para rodar com outro ambiente conda) por causa do inferno que é instalar o pacote 'terra' 
rule extract_variables:
    input: 
        meta = config['meta_path'],
        rasters = expand("../resources/climate_data/wc2.1_2.5m_{variable}_{year}-{month}.tif",  year = ["1997","2009","2010","2017"], month = ["{:02}".format(i) for i in range(1, 13)], variable = ["tmin", "tmax", "prec"])
    output: "../resources/climate_data/summary_{year}_6_months_prior.tsv"
    params: path = "../resources/climate_data"
    wildcard_constraints:  year = "(1997|2009|2010|2017)"
    conda: "/home/vitoria/time_clines_workflow/config/terra_env.yaml"
    shell: "Rscript ./scripts/R/organize_worldclim_data.R -meta {input.meta} -year {wildcards.year} -path {params.path} -o {output}"