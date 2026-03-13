#### MAGMA preprocess ####
require(future);require(furrr);require(purrr)

####### define kbp windows ########
pathology = c("SCZ","MDD","BIP","ASD","ADHD","OCD","PTSD","CUD","SA","AD","PD","ALS","CD","UC")
windows = c("35.10", "100")
target.dir = "script.data.enrichment/RData/"
output.dir = "McGill/seed_400/output/MAGMA/"
######## Annotation step ####
purrr::walk(windows,~{
  kbp = .x
  magma = "magma/magma"
  snp_loc  = paste0(target.dir,"g1000_eur.bim") ### snpmap of 1000 genome
  gene_loc = paste0(output.dir,"gene_loc") ### geneMap created previously
  if(kbp == "35.10"){window = paste0("window=",gsub("\\.",",",kbp))} else window = paste0("window=",kbp)
  out = paste0(output.dir,kbp,"kbp_annotation.step")
  cmd = paste0(magma," --annotate ",window," --snp-loc ",snp_loc," --gene-loc ",gene_loc," --out ",out)
  
  if (!grepl(paste0(kbp,"kbp\\_annotation\\.step\\.genes\\.annot"),list.files(output.dir))){
    ##### run analysis
    system(cmd)
  } else print(paste0(out," already present"))
})

####### Gene analysis step ####
purrr::walk(windows,~{
  kbp = .x
  purrr::walk(pathology,function(patho){
    sumstats = grep("check|clumped|grch38|\\.version|\\.RData|PGC2|PGC2\\.clozuk",
                    list.files(paste0(target.dir,tolower(patho),"/")),value = T,invert = T)   
    purrr::walk(sumstats,~{
      pgc = .x
      ### run Magma
      magma      = "magma/magma"
      gene_annot = paste0(output.dir,kbp,"kbp_annotation.step.genes.annot") ### output of annotation step
      pval       = paste0(target.dir,tolower(patho),"/",pgc," ncol=N") ### summary statistic of PGC
      out        = paste0(output.dir,kbp,"kbp_",pgc,"_gene.analysis")
      
      future::plan("multicore")
      dir.file = furrr::future_map_chr(1:22,~{
        chr = .x
        
        bfile = paste0(target.dir,"g1000_eur.chr",chr) ### 1000 genome genotype file for LD estimation
        cmd = paste0(magma," --bfile ",bfile," --gene-annot ",gene_annot," --pval ",pval," --batch ",chr," chr --out ",out)
        
        if ((patho == "RA" & chr == 6)|(patho == "UC" & chr == 6)) {
          cmd = paste0(cmd, " --gene-settings adap-permp=10000,10")
        }
        
        
        ##### run analysis
        raw.file = paste0(out, ".batch", chr, "_chr.genes.raw")
        # out.file = paste0(out, ".batch", chr, "_chr.genes.out")
        
        system(cmd)
        return(gsub(".genes.raw$","",raw.file))
        
      })
      future::plan("sequential")
      
      ### merge MAGMA batch output
      cmd.1 = paste0(paste0("{ cat ",dir.file[1],".genes.raw;"),paste0(" sed '1,2d' ",dir.file[-1],collapse = ".genes.raw;"),".genes.raw; } > ",
                     out,".genes.raw")
      cmd.2 = paste0(paste0("{ cat ",dir.file[1],".genes.out;"),paste0(" sed '1d' ",dir.file[-1],collapse = ".genes.out;"),".genes.out; } > ",
                     out,".genes.out")
      
      system(cmd.1)
      print(paste0("Writing file: ",out,".genes.raw"))
      system(cmd.2)
      print(paste0("Writing file: ",out,".genes.out"))
      
    })
  })
})

