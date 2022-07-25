
setwd('K:\\quercetin_AD_code_data')

set.seed(123456)

library(ggfun)
library(showtext)
library(sysfonts)
library(affy)
library(doMPI)
library(doParallel)
library(ggfortify)
library(biomaRt)
library(GEOquery)


eSet_GSE48350 = getGEO( 'GSE48350', destdir = './data',  # downloaded in [#02Annotation] folder
                       AnnotGPL = F, getGPL = F)  # remove platform info


eSet_GSE5281 = getGEO( 'GSE5281', destdir = './data',  # downloaded in [#02Annotation] folder
                       AnnotGPL = F, getGPL = F)  # remove platform info


save( eSet_GSE48350, file = './data/eSet_GSE48350.rdata')
save( eSet_GSE5281, file = './data/eSet_GSE5281.rdata')


load('./data/eSet_GSE48350.rdata')
load('./data/eSet_GSE5281.rdata')
pheno_GSE48350 = eSet_GSE48350[[1]]@phenoData@data
pheno_GSE5281 = eSet_GSE5281[[1]]@phenoData@data


pheno_GSE48350$batch_num = 1
pheno_GSE5281$batch_num = 2

library(plyr)
pheno = plyr::rbind.fill( pheno_GSE48350, pheno_GSE5281 )
View(pheno)



exp_GSE48350 = eSet_GSE48350[[1]]@assayData[["exprs"]]
exp_GSE5281 = eSet_GSE5281[[1]]@assayData[["exprs"]]


GPL_GSE48350 = eSet_GSE48350[[1]]@annotation
GPL_GSE5281 = eSet_GSE5281[[1]]@annotation


library(idmap2)
annotation_frame = get_soft_IDs('GPL570')
colnames(annotation_frame) = c('probe_id', 'symbol')




exp_GSE48350 = as.data.frame(exp_GSE48350)
exp_GSE48350$probe_id = rownames(exp_GSE48350)
exp_GSE48350$probe_id = gsub('(+)', '', exp_GSE48350$probe_id)
exp_GSE48350$probe_id = gsub(' ', '', exp_GSE48350$probe_id)
exp_GSE48350 = merge(exp_GSE48350, annotation_frame , by.x = 'probe_id', by.y = 'probe_id')
exp_GSE48350 = exp_GSE48350[, - which( colnames(exp_GSE48350)=='probe_id') ]
exp_GSE48350 = exp_GSE48350[which( !is.na(exp_GSE48350$symbol) ), ]

### parse multiple probes corresponding to the same gene
median_data = apply(exp_GSE48350[, colnames(exp_GSE48350)!='symbol'], 1, median)   
exp_GSE48350$median = median_data
exp_GSE48350 = exp_GSE48350[order(exp_GSE48350$symbol, exp_GSE48350$median, decreasing = T),]    
exp_GSE48350 = exp_GSE48350[!duplicated(exp_GSE48350$symbol),]   
rownames(exp_GSE48350) = exp_GSE48350$symbol
exp_GSE48350 = exp_GSE48350[, ! colnames(exp_GSE48350) %in% c('symbol', 'median') ]

anyNA(exp_GSE48350)
qx <- as.numeric(quantile(exp_GSE48350, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
isLog2_needed <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if(isLog2_needed){
    exp_GSE48350 = log2(exp_GSE48350+1)
}

save(exp_GSE48350, file='./data/exp_GSE48350_anotated.rdata')




exp_GSE5281 = as.data.frame(exp_GSE5281)
exp_GSE5281$probe_id = rownames(exp_GSE5281)
exp_GSE5281$probe_id = gsub('(+)', '', exp_GSE5281$probe_id)
exp_GSE5281$probe_id = gsub(' ', '', exp_GSE5281$probe_id)
exp_GSE5281 = merge(exp_GSE5281, annotation_frame , by.x = 'probe_id', by.y = 'probe_id')
exp_GSE5281 = exp_GSE5281[, - which( colnames(exp_GSE5281)=='probe_id') ]
exp_GSE5281 = exp_GSE5281[which( !is.na(exp_GSE5281$symbol) ), ]

### parse multiple probes corresponding to the same gene
median_data = apply(exp_GSE5281[, colnames(exp_GSE5281)!='symbol'], 1, median)   
exp_GSE5281$median = median_data
exp_GSE5281 = exp_GSE5281[order(exp_GSE5281$symbol, exp_GSE5281$median, decreasing = T),]    
exp_GSE5281 = exp_GSE5281[!duplicated(exp_GSE5281$symbol),]   
rownames(exp_GSE5281) = exp_GSE5281$symbol
exp_GSE5281 = exp_GSE5281[, ! colnames(exp_GSE5281) %in% c('symbol', 'median') ]

anyNA(exp_GSE5281)
qx <- as.numeric(quantile(exp_GSE5281, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
isLog2_needed <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if(isLog2_needed){
    exp_GSE5281 = log2(exp_GSE5281+1)
}

save(exp_GSE5281, file='./data/exp_GSE5281_anotated.rdata')







exp_GSE48350$symbol = rownames(exp_GSE48350)
exp_GSE5281$symbol = rownames(exp_GSE5281)
exp_GSE48350_GSE5281_merge = merge( exp_GSE48350, exp_GSE5281 , by.x = 'symbol', by.y = 'symbol')
rownames(exp_GSE48350_GSE5281_merge) = exp_GSE48350_GSE5281_merge$symbol
exp_GSE48350_GSE5281_merge = exp_GSE48350_GSE5281_merge[, - which(colnames(exp_GSE48350_GSE5281_merge)=='symbol')]

save( exp_GSE48350_GSE5281_merge, file='./data/exp_GSE48350_GSE5281_merge.rdata')

exp = exp_GSE48350_GSE5281_merge



reg_str_extract = function(str, reg){
    library(stringr)
    new_str = str_extract(str, reg)
    return(new_str)
}

one_colname = function(My_data_frame, a_colname, chooseit=1){
    if(  a_colname %in% colnames(My_data_frame)  ){
        return(    My_data_frame[, which(colnames(My_data_frame) %in% a_colname )]    )
    }else{
        return(My_data_frame[,grep(a_colname, colnames(My_data_frame))[chooseit]])
    }
}


get_group = function(  pheno, group_by_word='source_name', reg='', CaseName='', ControlName='', reg2Case='', reg2Control='',  group_no_numeric=F, suffix_no_numeric=T, is_replace_space=F, add_columns=''   ){
    gorup_frame = transform(pheno, sample = rownames(pheno))
    pheno_temp = pheno
    ControlName = gsub(' ', '_', ControlName)
    CaseName = gsub(' ', '_', CaseName)
    
    if( reg=='' & reg2Case=='' & reg2Control=='' ){
        if(group_no_numeric==T){
            group_temp = data.frame(rownames(pheno_temp), gsub("\\d","",one_colname(pheno_temp, group_by_word)), pheno_temp$batch_num )}
        else{
            if(suffix_no_numeric==T){
                group_temp = data.frame(rownames(pheno_temp), sub('[0-9]+$', '', one_colname(pheno_temp, group_by_word)), pheno_temp$batch_num )
            }else{
                group_temp = data.frame(rownames(pheno_temp), one_colname(pheno_temp, group_by_word), pheno_temp$batch_num )
            }
        }
    }
    else{
        if( reg!='' & reg2Case=='' & reg2Control==''  ){
            group_temp = data.frame( rownames(pheno_temp),  reg_str_extract( one_colname(pheno_temp, group_by_word), reg ),  pheno_temp$batch_num )
            colnames(group_temp) = c('sample', 'group', 'batch_num')
        }else if(   reg2Case!='' & reg2Control!='' & reg==''  ){
            reg = paste0( reg2Case, '|', reg2Control )
            reg = gsub('\\|\\|', '', reg)
            group_temp = data.frame( rownames(pheno_temp),  reg_str_extract( one_colname(pheno_temp, group_by_word), reg ),  pheno_temp$batch_num )
            colnames(group_temp) = c('sample', 'group', 'batch_num')
            if( CaseName != '' ){
                group_temp$group = sub( reg2Case, CaseName, group_temp$group )
            }else{
                group_temp$group = sub( reg2Case, "Case", group_temp$group )
            }
            if( ControlName != '' ){
                group_temp$group = sub( reg2Control, ControlName, group_temp$group )
            }else{
                group_temp$group = sub( reg2Control, "Control", group_temp$group )
            }
            
        }else{
            print('-=-=-=-【请输入正确的reg、reg2Case、reg2Control】-=-=-=-')
        }
        
    }
    
    colnames(group_temp) = c('sample', 'group', 'batch_num')
    if(!is.data.frame(group_temp)){group_temp = as.data.frame(group_temp)}
    rownames(group_temp) = group_temp$sample
    
    if(is_replace_space){
        group_temp$group = gsub(' ', '_', group_temp$group)
    }
    
    if( add_columns!='' ){
        for(a_col in add_columns){
            if(a_col != ''){
                group_temp[, a_col] = pheno[rownames(pheno) %in% group_temp$sample, a_col] 
            }
        }
    }
    return( group_temp )
}



pheno$bath_num = pheno$batch_num
rownames(pheno) = pheno$geo_accession
pheno$group_by = paste0(pheno$title, '@', pheno$source_name_ch1)
# group = get_group(  pheno, group_by_word='group_by', reg='', CaseName='', ControlName='', reg2Case="(^Frontal_Cortex AD)|(^superior[ ]*frontal[ ]*gyrus_[a-zA-Z]+_[0-9]+_AD_)|(^SFG_affected)|(^Subject[ 0-9], region P*r*e*[fF]*rontal)|(^AD_FC[^@]+@Frontal cortex)", reg2Control='(^Frontal_Cortex control)|(^SuperiorFrontalGyrus_[a-zA-Z]+_[0-9]+yrs_indiv)|(^SFG control)|(^non-AD_FC[^@]+@Frontal cortex)', group_no_numeric=F, suffix_no_numeric=T   )#额叶
group = get_group(  pheno, group_by_word='group_by', reg='', CaseName='', ControlName='', reg2Case="(^Frontal_Cortex AD)|(^superior[ ]*frontal[ ]*gyrus_[a-zA-Z]+_[0-9]+_AD_)|(^SFG_affected)|(^Subject[ 0-9], region P*r*e*[fF]*rontal)|(^AD_FC[^@]+@Frontal cortex)", reg2Control='(^Frontal_Cortex control)|(^SuperiorFrontalGyrus_[a-zA-Z]+_[0-9]+yrs_indiv)|(^SFG control)|(^non-AD_FC[^@]+@Frontal cortex)', group_no_numeric=F, suffix_no_numeric=T   )#额叶
# group = get_group(  pheno, group_by_word='title', reg='', CaseName='',ControlName='', reg2Case="(_affected_)|(_AD_)", reg2Control='( control )|(yrs_indiv)', group_no_numeric=F, suffix_no_numeric=T   )
# colnames(group)[colnames(group)=='bath_num'] = 'batch_num'






exp_in_group = function(  exp, group, chose_list='', order_col='sample', is_decreasing=F  ){
    if(!is.data.frame(exp)){ exp = as.data.frame( exp ) }
    if(!is.data.frame(group)){ group = as.data.frame( group ) }
    
    
    # colnames(group) = c('sample', 'group', colnames(group)[3])
    
    intersection_list1 = intersect(colnames(exp), rownames(group))
    # intersection_list2 = intersect(intersection_list1, chose_list)
    
    
    if(  chose_list!= '' ){
        group_chosed = group[which(group$group %in% chose_list), ] 
    }else{
        group_chosed = group
    }
    group_chosed = group_chosed[which(group_chosed$sample %in% intersection_list1), ]
    
    
    if(  order_col!= ''   &   order_col %in%  colnames(group_chosed) ){
        print(paste0('正在按[', order_col, ']排序！！！'))
        group_chosed = group_chosed[ order(group_chosed[, order_col], decreasing=is_decreasing) ,]
    }
    
    
    colnames(exp) = gsub('\\.x', '', colnames(exp))
    ### 2个相交
    # exp_chosed = exp[, which(colnames(exp) %in% group_chosed$sample )]
    ### 
    exp_chosed = exp[,   colnames(exp)[match(  group_chosed$sample , colnames(exp)  )]          ]
    
    
    # exp_chosed = exp[, which(colnames(exp) %in% chose_list )]
    # sample_list = NULL
    # for(i in colnames(exp)){
    #   sample_list = c(sample_list, i %in% group$sample)
    # }
    # exp = exp[, sample_list]
    
    # group_chosed = group_chosed[ order(group_chosed$sample)    ,]
    # exp_chosed = exp_chosed[, order(colnames(exp_chosed))]
    return(   list(exp_chosed, group_chosed)     )
}



not_NA_group_array = unique(group$group)[!is.na(unique(group$group))]
exp_in_group_list = exp_in_group(   exp, group, chose_list = not_NA_group_array     )
exp = exp_in_group_list[[1]]
group = exp_in_group_list[[2]]





library(sva)
exp_removed_batch_effect = data.frame(   ComBat(exp, group$batch_num, model.matrix(~as.factor(   group$group   )), par.prior=TRUE)    )

save( exp_removed_batch_effect, file='./data/exp_removed_batch_effect.rdata')




# load( './data/exp_removed_batch_effect.rdata' )
boxplot( exp_removed_batch_effect, boxwex=0.7, notch=T, main='test', outline=FALSE, las=2)

Pett_PCA_plot = function( exp, group_frame, GSE_or_Sample='Sample', method=0 ){
    exp = as.data.frame( exp )
    
    
    if(   length(table(group_frame$group))  > 20    ){
        print(     paste0('---【group$group分组数量太多，[ >20 ]，提前终止！】---')      )
        return(   paste0('---【group$group分组数量太多，[ >20 ]，提前终止！】---')   )
    }
    
    # GSE_or_Sample='Sample' ### 或者 GSE_or_Sample='GSE'
    method = as.numeric( method )
    if(F){library(ggpubr);show_point_shapes()}  ###### 显示点的形状(ggpubr包内的函数）
    
    library(ggplot2)
    if( dim(exp)[1]==dim(group_frame)[1] ){
        exp_t = as.data.frame( exp )
    }else{
        print(' t() ')
        exp_t = as.data.frame( t(exp) )
    }
    dim(exp_t)
    
    dim(group_frame)
    Batch_num = as.character(group_frame$batch_num)
    
    if( method!=2 ){
        library("FactoMineR")   ##### 方法 1
        library("factoextra")   ##### 方法 1
        pca_FactoMineR <- PCA( exp_t, graph = FALSE)   ##### 方法 1
        if(GSE_or_Sample=='Sample'){ pca_FactoMineR_res = get_pca_ind( pca_FactoMineR ) }else if(GSE_or_Sample=='GSE'){ pca_FactoMineR_res = get_pca_var( pca_FactoMineR ) }else{ print('---【[GSE_or_Sample]参数错误】---') }
        pca_original_data = pca_FactoMineR_res$coord
        pca_original_data = as.data.frame(pca_original_data)
    }
    else{
        pca_res1 <- prcomp( exp_t )    ##### 方法 2
        if(GSE_or_Sample=='Sample'){ pca_original_data = data.frame(pca_res1$x) }else if(GSE_or_Sample=='GSE'){ pca_original_data = data.frame(pca_res1$rotation) }else{ print('---【[GSE_or_Sample]参数错误】---') }
    }
    
    colnames(pca_original_data)[1:2] = c('PC1', 'PC2')
    
    # pcadata = pca_res1$rotation
    # pcadata<-data.frame(sample=rownames(pcadata),pcadata) %>% left_join(group,by="sample")
    group = factor(group_frame$group)
    batch = factor(Batch_num)
    png = ggplot(pca_original_data, aes(x=PC1, y=PC2)) + 
        # geom_hline(yintercept = 0,linetype="dashed") + # 添加横线
        # geom_vline(xintercept = 0,linetype="dashed") + # 添加竖线
        theme_bw() + # 加上边框
        # stat_ellipse(aes(color = Group), level = 0.95, linetype = 2, show.legend = T) + # 添加置信椭圆
        geom_point(aes(shape=batch, color=group), size=3) +  # 画散点图并设置大小
        scale_shape_manual(values = c(3, 1, 17, 16, 2, 10, 11, 6, 7, 0, 9, 5, 12, 13, 14, 15,  18, 19, 8, 4))  #自定义点的形状，分别为15， 19， 17
    # stat_ellipse(level=0.95,show.legend=F,geom="polygon",alpha=0.2)
    
    png
    
    return(list(pca_original_data, png))
}

pca_data1 = Pett_PCA_plot( exp_removed_batch_effect, group, GSE_or_Sample='Sample', method=0 )
pca_original_data1 = pca_data1[[1]]
png_pca1 = pca_data1[[2]]
png_pca1



exp_normalize_quantiles = function(exp){
    library("preprocessCore")
    exp_norm = normalize.quantiles(as.matrix(exp), copy=TRUE)
    # boxplot(exp_norm, boxwex=0.7, notch=T, main=MyGSE, outline=FALSE, las=2)
    colnames(exp_norm) = colnames(exp)
    rownames(exp_norm) = rownames(exp)
    # tryCatch({    rm(c('exp', 'exp_norm'))           }, error=function(cond){        11111        })
    # boxplot(exp, boxwex=0.7, notch=T, main=MyGSE, outline=FALSE, las=2)
    return(as.data.frame(exp_norm))
}

exp_quantiles = exp_normalize_quantiles( exp_removed_batch_effect ) ###############################################
boxplot( exp_quantiles, boxwex=0.7, notch=T, main='test', outline=FALSE, las=2)


pca_data2 = Pett_PCA_plot( exp_quantiles, group, GSE_or_Sample='Sample', method=0 )
pca_original_data2 = pca_data2[[1]]
png_pca2 = pca_data2[[2]]
png_pca2








not_NA_group_array = unique(group$group)[!is.na(unique(group$group))]
exp_in_group_list = exp_in_group(   exp_quantiles, group, chose_list = not_NA_group_array     )
exp = exp_in_group_list[[1]]
group = exp_in_group_list[[2]]




library(limma)
library(doMPI)
library(doParallel)
library(DESeq2)
library(edgeR)




control = 'Control'
##### 进行差异分析
library(limma)
# 2.做分组矩阵
design <- model.matrix(~0+factor(   one_colname(group, 'group')   ))
colnames(design) = levels(factor(   one_colname(group, 'group')   ))
# rownames(design)=colnames(exp)
head(design)      #分组矩阵

# 3.做比较矩阵
if( control!='' ){
    if(   colnames(design)[2]!= control  ){
        temp_array = design[, 2]
        design[, 2] = design[, 1]
        design[, 1] = temp_array
    }
}
colnames(design) = c('Case', 'Control')
contrast.matrix = makeContrasts(paste0(c('Case', 'Control'), collapse = "-"),levels = design)
# contrast.matrix<-makeContrasts(paste0(c(colnames(design)[1], colnames(design)[2]),collapse = "-"),levels = design)
contrast.matrix   #差异比较矩阵
#到此，做差异分析所需要的三个矩阵就做好了：表达矩阵、分组矩阵、差异比较矩阵#
#我们已经制作好了必要的输入数据，下面开始讲如何使用limma这个包来进行差异分析了！

##step1
fit <- lmFit(exp, design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果#
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
#找出差异基因检验结果并输出符合条件的结果  这里更改pvalue和logFC
DEGs_df = topTable(fit2, adjust="fdr", coef=1, n=Inf, sort.by = 'logFC')


ADGs_df = subset(ADGs_df,   abs(logFC) > log2(1.3)   &   adj.P.Val < 0.05    )
save( ADGs_df, file='./data/ADGs_df.rdata')

library('xlsx')
xlsx::write.xlsx2( ADGs_df, file='./data/ADGs_df.xlsx', sheetName='ADGs', showNA = TRUE )

### upregulated genes
dim( subset(ADGs_df,   logFC > log2(1.3)   &   adj.P.Val < 0.05    ))
### downregulated genes
dim( subset(ADGs_df,   logFC < -log2(1.3)   &   adj.P.Val < 0.05    ))





# ADGs = Pett_read_write_table("E:\\_____writing_____\\002_drug_immune_target_Quercetin\\补充材料原来\\Supplementary Material Table 2.xlsx")
# ADGs = as.data.frame(ADGs)
# rownames(ADGs) = ADGs$gene_symbol
# ADGs = ADGs[, - which(colnames(ADGs)=='gene_symbol')]


# save.image('暂时002.rdata')




load('./data/ADGs_df.rdata')
library(DOSE)
library(clusterProfiler)#进行GO富集和KEGG富集
library('org.Hs.eg.db') #人类注释数据库
library('org.Mm.eg.db') #小鼠注释数据库
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(enrichplot)# 可以绘制【GSEA】
library(ReactomePA)
library(doMPI)




id_convert_symbol_ENTREZID = function( DEG_frame, Taxonomy='homo',  symbol_col_name='probe_id', fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL","SYMBOL", "UNIPROT")   ){
    ### 'fromType' should be one of ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, 
    ### ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, 
    ### GENETYPE, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, 
    ### PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIPROT.
    
    fromType = toupper(fromType)
    toType = toupper(toType)
    library(clusterProfiler)#进行GO富集和KEGG富集
    library('org.Hs.eg.db') #人类注释数据库
    library('org.Mm.eg.db') #小鼠注释数据库
    # Taxonomy目前只有人类和小鼠
    if(   class(DEG_frame) == 'character'   ){
        print('--【输入是array】--')
        DEG_frame = data.frame(gene_input = DEG_frame)
        symbol_col_name='gene_input'
    }else{
        print('--【输入是data.frame】--')
        DEG_frame = data.frame( DEG_frame )
    }
    gene_input = gene_array_get_duplicated_gene( DEG_frame[symbol_col_name][,1], save_old_duplicated_gene_name=F)
    OrgDb = ifelse(tolower(Taxonomy) %in% c('homo', 'homo sapiens', 'human', '9606', 9606), "org.Hs.eg.db", 'org.Mm.eg.db')
    symbol_ENTREZID = bitr( gene_input, fromType=fromType,
                            toType=unique(c(fromType, toType)),
                            OrgDb=OrgDb) # 24.65% of input gene IDs are fail to map...
    DEG_frame = left_join_dt_Pett(DEG_frame, symbol_ENTREZID, by1=symbol_col_name, by2=fromType)
    return(DEG_frame)
}


gene_array_get_duplicated_gene = function(gene_array, save_old_duplicated_gene_name=T){
    duplicated_gene_name_array = gene_array[      grep('///', gene_array)       ]
    if(   length(duplicated_gene_name_array) > 0   ){
        
        
        all_duplicated_gene_array = NULL
        all_list = list()
        total_num = length(duplicated_gene_name_array)
        for(num in 1:total_num){
            Pett_ProgressBar(num, total_num)
            duplicated_gene_combined = duplicated_gene_name_array[num]
            for(a_name in strsplit(duplicated_gene_combined,split='///')[[1]]){
                if((! a_name %in% all_duplicated_gene_array)  &  (! a_name %in% rownames(exp))  ){
                    all_duplicated_gene_array = c(all_duplicated_gene_array, a_name)
                    all_list[[a_name]] = duplicated_gene_combined
                }
            }
        }
        
        
        gene_array_plus_duplicated_gene = c(gene_array, all_duplicated_gene_array)
        
        if( !save_old_duplicated_gene_name ){
            gene_array_plus_duplicated_gene = gene_array_plus_duplicated_gene[  -  grep('///', gene_array_plus_duplicated_gene)    ]
        }
        
        return(gene_array_plus_duplicated_gene)
        
    }else{
        print('        --【没有含有///的基因】--')
        return(   gene_array   )
    }
    
}


Pett_ProgressBar = function(start_num, total_num, close=F ){
    style=3
    setTxtProgressBar( txtProgressBar(style=style), start_num/total_num )
    # close( txtProgressBar(style=style) )
    # 
    # if( start_num < (total_num) ){
    #   setTxtProgressBar( txtProgressBar(style=style), start_num/total_num )
    # }else{
    #   setTxtProgressBar( txtProgressBar(style=style), 1 )
    #   setTxtProgressBar( txtProgressBar(style=style), start_num/(total_num) )
    #   if( close ){ close( txtProgressBar(style=3) ) }
    # }
}

left_join_dt_Pett = function(table_or_path1, table_or_path2, by1, by2, colname2_prefix='', colname2_suffix='', remove_dulpulicate_colname_of_2='auto' ){
    library(tidyfst)
    library(data.table)
    # source('C:\\Users\\Administrator\\Documents\\#File#\\OneDrive\\Win10\\Python\\_____001_R_Modules_____\\001_GEO_module.R', encoding = 'UTF-8')
    table_or_path1_class = 'data.frame'
    if(length(class(table_or_path1_class))>1){table_or_path1_class='NOT_data.frame'}
    
    if(   (! 'data.frame' %in% class(table_or_path1))  &  class(table_or_path1) == 'character'  ){
        if(    grepl('`.`', table_or_path1)    ){
            table_or_path1 = gsub('^`', '', table_or_path1)
            table_or_path1 = gsub('`$', '', table_or_path1)
            table_or_path1 = strsplit(table_or_path1, split = "`.`")
            table_or_path1 = dbsql(  paste0("select  *   from `", table_or_path1[[1]][2], "`;"), table_name='', sql_or_litePath='sql', dbname = table_or_path1[[1]][1], limited_num=Inf,  append=T,  dbListTables_only=F  )
        }else{
            table_or_path1 = as.data.table(Pett_read_write_table(   table_or_path1   )) 
        }
    }
    
    if(   (! 'data.frame' %in% class(table_or_path2))  &  class(table_or_path2) == 'character'  ){
        if(    grepl('`.`', table_or_path2)    ){
            table_or_path2 = gsub('^`', '', table_or_path2)
            table_or_path2 = gsub('`$', '', table_or_path2)
            table_or_path2 = strsplit(table_or_path2, split = "`.`")
            table_or_path2 = dbsql(  paste0("select  *   from `", table_or_path2[[1]][2], "`;"), table_name='', sql_or_litePath='sql', dbname = table_or_path2[[1]][1], limited_num=Inf,  append=T,  dbListTables_only=F  )
        }else{
            table_or_path2 = as.data.table(Pett_read_write_table(   table_or_path2   )) 
        }
    }
    
    
    if(tolower(by1) %in% c('rowname', 'rownames', '行', '行名')){
        table_or_path1$Pett_temp_by1_column_12345678 = rownames(table_or_path1)
        by1 = 'Pett_temp_by1_column_12345678'
    }
    if(tolower(by2) %in% c('rowname', 'rownames', '行', '行名')){
        table_or_path2$Pett_temp_by2_column_12345678 = rownames(table_or_path2)
        by2 = 'Pett_temp_by2_column_12345678'
    }
    # if(   ! 'data.frame' %in% class(table_or_path2)   ){
    #   table_or_path2 = as.data.table(Pett_read_write_table(   table_or_path2   ))
    # }
    if(   ! 'data.table' %in% class(table_or_path1)   ){
        table_or_path1 = as.data.table(table_or_path1)
    }
    if(   ! 'data.table' %in% class(table_or_path2)   ){
        table_or_path2 = as.data.table(table_or_path2)
    }
    
    
    ### ### ### ### ### ### ### ### ###
    if( length(by1) > 1 ){
        # temp_combined_colname = paste0(by1, collapse='+++')
        # temp_combined_Value = NULL
        for(num in 1:length(by1)){
            setnames(table_or_path1, c(by1[num]), c('temp_test_colname_of_by1_1008611'))
            if(num==1){     
                table_or_path1 =  table_or_path1[, temp_test_colname_of_by1_total := temp_test_colname_of_by1_1008611    ]
                setnames(table_or_path1, c('temp_test_colname_of_by1_1008611'), c(by1[num]))
            }else{
                table_or_path1 =  table_or_path1[, temp_test_colname_of_by1_total := paste0(temp_test_colname_of_by1_total, '+++', temp_test_colname_of_by1_1008611)    ]     
                setnames(table_or_path1, c('temp_test_colname_of_by1_1008611'), c(by1[num]))
            }
        }
        by1 = 'temp_test_colname_of_by1_total'
    }
    
    ### ### ### ### ### ### ### ### ###
    if( length(by2) > 1 ){
        for(num in 1:length(by2)){
            setnames(table_or_path2, c(by2[num]), c('temp_test_colname_of_by2_1008611'))
            if(num==1){     
                table_or_path2 =  table_or_path2[, temp_test_colname_of_by2_total := temp_test_colname_of_by2_1008611    ]
                table_or_path2[, temp_test_colname_of_by2_1008611 := NULL ]
            }else{   
                table_or_path2 =  table_or_path2[, temp_test_colname_of_by2_total := paste0(temp_test_colname_of_by2_total, '+++', temp_test_colname_of_by2_1008611)    ]     
                table_or_path2[, temp_test_colname_of_by2_1008611 := NULL ]
            }
        }
        by2 = 'temp_test_colname_of_by2_total'
    }
    
    
    ### 添加【前缀】、【后缀】
    real_colname_of_table_or_path2 = by2
    if( colname2_prefix != ''){
        colnames(table_or_path2) = paste0(colname2_prefix , colnames(table_or_path2))
        real_colname_of_table_or_path2 = paste0(colname2_prefix , real_colname_of_table_or_path2)
    }
    if( colname2_suffix != ''){
        colnames(table_or_path2) = paste0(  colnames(table_or_path2), colname2_suffix  )
        real_colname_of_table_or_path2 = paste0( real_colname_of_table_or_path2, colname2_suffix  )
    }
    colnames(table_or_path2)[colnames(table_or_path2)==real_colname_of_table_or_path2] = by1
    
    
    ### ### ### ### ### ### ### ### ### 【[合并]】 ### ### ### ### ### ### ### ### ###
    table_or_path1 = left_join_dt(table_or_path1, table_or_path2, by1)
    
    ### 去掉多余的列
    if(   'temp_test_colname_of_by1_total' %in%  colnames(table_or_path1)   ){    table_or_path1[, temp_test_colname_of_by1_total := NULL ]               }
    if(   'Pett_temp_by1_column_12345678' %in%  colnames(table_or_path1)   ){    rownames(table_or_path1)=table_or_path1$Pett_temp_by1_column_12345678;table_or_path1[, Pett_temp_by1_column_12345678 := NULL ]               }
    if(   'Pett_temp_by2_column_12345678' %in%  colnames(table_or_path1)   ){    
        # rownames(table_or_path1)=table_or_path1$Pett_temp_by2_column_12345678;
        table_or_path1[, Pett_temp_by2_column_12345678 := NULL ]               }
    # table_or_path1$Pett_temp_by1_column_12345678
    if(table_or_path1_class == 'data.frame'){ table_or_path1$last_one_temp_123456 = rownames(table_or_path1);table_or_path1 = as.data.frame(table_or_path1); rownames(table_or_path1)=table_or_path1$last_one_temp_123456;table_or_path1=table_or_path1[, -dim(table_or_path1)[2]]   }
    return(table_or_path1)
}
convert_symbol = symbol_convert = gene_convert_symbol_ENTREZID = id_convert_symbol_ENTREZID


DEG_symbol_ENTREZID = gene_convert_symbol_ENTREZID(unique(rownames(ADGs_df)))

DEGs_chosed_ENTREZID_array = DEG_symbol_ENTREZID$ENTREZID


library('org.Hs.eg.db') #人类注释数据库
library('org.Mm.eg.db') #小鼠注释数据库
BP_MF_CC_ALL = 'BP'
BP_MF_CC_ALL = 'ALL'
Taxonomy = 'homo'
DEG_FDR_cutoff = 0.05
GO_p.val_cutoff = 1
GO_q.val_cutoff = 1

OrgDb = ifelse(tolower(Taxonomy) %in% c('homo', 'homo sapiens', 'human', '9606', 9606), "org.Hs.eg.db", 'org.Mm.eg.db')
ego <- enrichGO(
    # gene = DEGs_chosed_ENTREZID_array[!is.na(DEGs_chosed_ENTREZID_array)],
    # gene = all_gene_ENTREZID_array,
    # gene = all_input_gene_ENTREZID_array ,
    gene = DEGs_chosed_ENTREZID_array ,
    # universe      = DEG_symbol_ENTREZID$ENTREZID,
    OrgDb         = OrgDb,
    ont           = toupper(BP_MF_CC_ALL), ### One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    # ont           = toupper('ALL'), ### One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    ### BP（生物学过程）、MF（分子功能）和CC（细胞组成）
    ### 
    ### 分子功能(Molecular Function,MF )
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_p.val_cutoff,
    qvalueCutoff  = GO_q.val_cutoff,
    readable      = TRUE)
ego_result_DF = ego@result

png_GO_barplot = barplot(ego, title = "Biological process", showCategory=20)

library('cowplot')
png(file='Figure 3A.png', height=600, width=400, bg="white")
print(png_GO_barplot)
dev.off()





get_list_names = function(  this_list   ){
    return(      names(sapply(this_list, names))      )
}
list_key_array_get = get_list_names

get_date = function(format_date="%Y-%m-%d"){
    return(     format(Sys.time(), format=format_date)     )
}
get_time = function(format_date="%H:%M:%S", for_write_file_in_windows=F){
    if(for_write_file_in_windows){
        return(     format(Sys.time(), format="%H=%M=%S")     )
    }else{
        return(     format(Sys.time(), format=format_date)     )
    }
}
plot_GOChord = function( ego, DEG_frame,  GO_KEGG_chosed_Term_array = NA, max_num_Term_num = 10, max_num_gene_num = 50,      space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5 , png_height=1200, png_width=1100, png_backgroud="transparent"    ){
    library(GOplot)
    library(stringr)
    
    if(   is.na(GO_KEGG_chosed_Term_array) | GO_KEGG_chosed_Term_array == ''    ){
        # print('is.na(GO_KEGG_chosed_array) → NA')
        chord_prepare_DF = ego@result
        chord_prepare_DF = chord_prepare_DF[  1:max_num_Term_num   ,]
    }else{
        chord_prepare_DF = ego@result
        chord_prepare_DF = chord_prepare_DF[  chord_prepare_DF$Description %in% GO_KEGG_chosed_Term_array   ,]
    }
    
    all_term_array = chord_prepare_DF$Description
    
    
    GO_KEGG_chosed_Term_list = list()
    total_num = dim(chord_prepare_DF)[1]
    for(num in 1:total_num){
        # if(      num > max_num_Term_num    ){break;}
        
        Pett_ProgressBar(num, total_num)
        
        this_Term = chord_prepare_DF[   num   ,]$Description
        
        gene_array0 = chord_prepare_DF[   num   ,]$geneID
        English_len = nchar(   gsub('[0-9/]+', '', gene_array0)   )
        gene_array = strsplit(gene_array0, '/')[[1]]
        if(    English_len < length(gene_array)     ){
            ### 把【ENTREZID】转换成【SYMBOL】
            print('---【English_len <= 3】，把【ENTREZID】转换成【SYMBOL】---')
            gene_array_new = gene_convert_symbol_ENTREZID( data.frame( this_ENTREZID = gene_array ), Taxonomy='homo',  symbol_col_name='this_ENTREZID', fromType="ENTREZID", toType=c("ENTREZID","ENSEMBL","SYMBOL", "UNIPROT")   )$SYMBOL
            if(length(gene_array_new)<1){
                gene_array_new = gene_convert_symbol_ENTREZID( data.frame( this_ENTREZID = gene_array ), Taxonomy='mus',  symbol_col_name='this_ENTREZID', fromType="ENTREZID", toType=c("ENTREZID","ENSEMBL","SYMBOL", "UNIPROT")   )$SYMBOL
            }
            gene_array = gene_array_new
        }
        
        for(this_gene in gene_array){
            GO_KEGG_chosed_Term_list[[  this_gene   ]] = c( GO_KEGG_chosed_Term_list[[  this_gene   ]] ,  this_Term  )
        }
    }
    
    
    # test = data.frame( GO_KEGG_chosed_Term_list )
    # test = data.frame( GO_KEGG_chosed_Term_list )
    all_gene_array = get_list_names(GO_KEGG_chosed_Term_list)
    total_num = length(all_gene_array)
    GO_KEGG_chosed_gene_list = list()
    for(num in 1:total_num){
        this_gene = all_gene_array[num]
        this_Term_array = GO_KEGG_chosed_Term_list[[this_gene]]
        
        defult_value_array = rep(0, length(all_term_array))
        defult_value_array[which(all_term_array %in% this_Term_array)] = 1
        
        GO_KEGG_chosed_gene_list[[ this_gene ]] = defult_value_array
    }
    
    
    gene_Term_DF = as.data.frame(t(data.frame(GO_KEGG_chosed_gene_list)))
    colnames(gene_Term_DF) = all_term_array
    gene_Term_DF$symbol000 = rownames(gene_Term_DF)
    
    DEG_frame$probe_id = rownames(DEG_frame)
    DEG_frame = as.data.frame(DEG_frame)
    gene_Term_DF = as.data.frame(left_join_dt_Pett(gene_Term_DF, DEG_frame, by1='symbol000', by2='probe_id'))
    rownames(gene_Term_DF) = gene_Term_DF$symbol000
    
    chord_DF = gene_Term_DF[, colnames(gene_Term_DF) %in% c(all_term_array, 'logFC')]
    if(  dim(chord_DF)[1]  > max_num_gene_num    ){
        chord_DF = chord_DF[   1:max_num_gene_num   ,]
    }
    
    
    
    library(showtext) # 加载包
    library(sysfonts) # 加载包
    font_add("Times_Pett", "times.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Timesbd_Pett", "timesbd.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Timesi_Pett", "timesi.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Timesbi_Pett", "timesbi.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Arial_Pett", "arial.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Arialbd_Pett", "arialbd.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Ariali_Pett", "ariali.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_add("Arialbi_Pett", "arialbi.ttf") # myFont1赋予字体的名称，timesbd.ttf 为 Times New Roman粗体
    font_families() # 查看当前可用的字体名称
    file_name = paste0("GOChord[", get_date(), ' ', get_time(for_write_file_in_windows=T),"].png")
    png(file=file_name, height=png_height, width=png_width, bg=png_backgroud) #white, 
    showtext_begin()
    
    library(ggfun)
    # png_vocalno_set_font = set_font(png_vocalno, family="Helvetica", fontface="italic", color='firebrick')
    # png_vocalno_set_font = set_font(png_vocalno, family="myFont1", fontface='bold')
    
    # colnames(chord_DF)[colnames(chord_DF)=='logFC'] = 'log2(FoldChange)'
    # chord_DF$logFC = 2^chord_DF$logFC
    # colnames(chord_DF)[colnames(chord_DF)=='logFC'] = 'FoldChange'
    # tryCatch({      rm('ego', 'DEG_frame', 'gene_Term_DF', 'chord_DF', 'chord_prepare_DF')         }, error=function(cond){         111       })
    png_GOChord = GOChord( chord_DF, space = space, gene.order = 'logFC', gene.space = gene.space, gene.size = gene.size  )
    # png_GOChord = set_font(png_GOChord, family="Arialbd_Pett" )
    png_GOChord = set_font(png_GOChord, family="Arialbd_Pett" )
    png_GOChord
    # print('---【图片，请到工作目录找pdf文件】---')
    # print(paste0('~~文件名为【', file_name, '】~~'))
    ## 关闭showtext字体设置
    showtext_end() # 绘图命令放在showtext_begin()和showtext_end()之间
    # showtext_auto(FALSE) # 【[showtext_auto()]全局使用】不在需要就关闭
    dev.off()
    
    
    return(   png_GOChord   )
}



png_GOChord = plot_GOChord( ego, ADGs,  GO_KEGG_chosed_Term_array = NA, max_num_Term_num = 10, max_num_gene_num = 50,      space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5, png_height=1450, png_width=1400, png_backgroud='white'     )

library('cowplot')
png(file='Figure 3B.png', height=600, width=400, bg="white")
print(png_GO_barplot)
dev.off()






KEGG_p.val_cutoff = 1
KEGG_q.val_cutoff = 1
organism = ifelse(tolower(Taxonomy) %in% c('homo', 'homo sapiens', 'human', '9606', 9606), "hsa", 'mmu')
print('Make sure the Internet is connected!')
kegg_data <- enrichKEGG(
    # gene = DEGs_chosed_ENTREZID_array[!is.na(DEGs_chosed_ENTREZID_array)], ### gene = gene_up$ENTREZID,
    # gene = all_gene_ENTREZID_array,
    gene = DEGs_chosed_ENTREZID_array,
    # universe      = DEG_symbol_ENTREZID$ENTREZID,
    # gene = symbol_ENTREZID$ENTREZID, ### gene = gene_up$ENTREZID,
    organism = organism, # homo→hsa；mice→mmu；
    keyType = "kegg", # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
    # universe = gene_all,
    pvalueCutoff = KEGG_p.val_cutoff,
    qvalueCutoff = KEGG_q.val_cutoff) #P值或者Q值阈值可以根据富集的结果进行主观修改

png_KEGG_dotplot = enrichplot::dotplot( kegg_data, title = "KEGG pathways", showCategory = 20, )#气泡图
png(file='Figure 3C.png', height=600, width=400, bg="white")
print(png_KEGG_dotplot)
dev.off()



png_KeggChord = plot_GOChord( kegg_data, DEG_chosed_all,  GO_KEGG_chosed_Term_array = NA, max_num_Term_num = 10, max_num_gene_num = 50,      space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5, png_height=1450, png_width=1400, png_backgroud='white'     )
png(file='Figure 3D.png', height=600, width=400, bg="white")
print(png_KEGG_dotplot)
dev.off()




#####Immune Infiltration#####
library('CIBERSORT')
source("./code/CIBERSORT_Modules.R", encoding = 'UTF-8')
load('./data/LM22.rda')
TME.results = cibersort_run(LM22, exp, perm = 1000, QN = TRUE)

save(TME.results, file='./data/TME.results.rdata')




######################################################################################
######################################################################################
plot_vioplot_immune = function(TME.results, group, abscissa_distance=0.03 ){
    library('vioplot')
    
    rt = TME.results[, 1:22 ]
    tryCatch({      rm('TME.results')         }, error=function(cond){        111        })
    
    get_date = function(format_date="%Y-%m-%d"){
        return(     format(Sys.time(), format=format_date)     )
    }
    get_time = function(format_date="%H:%M:%S", for_write_file_in_windows=F){
        if(for_write_file_in_windows){
            return(     format(Sys.time(), format="%H=%M=%S")     )
        }else{
            return(     format(Sys.time(), format=format_date)     )
        }
    }
    
    
    file_name = paste0("vioplot_immune[", get_date(), ' ', get_time(for_write_file_in_windows=T),"].pdf")
    pdf(  file_name,  height=8, width=15)              #保存图片的文件名称
    par(las=1,mar=c(10,6,3,3))
    x=c(1:ncol(rt))
    y=c(1:ncol(rt))
    plot(x,y,
         xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
         main="",xlab="", ylab="Fraction",
         pch=21,
         col="white",
         xaxt="n")
    
    ### 找出control组、Case组的index
    all_group_array = unique(group$group)
    group_control_name = 'Control'
    if( ! group_control_name %in% group$group ){
        group_control_name = all_group_array[1]
    }
    group_case_name = 'Case'
    if( ! group_case_name %in% group$group ){
        group_case_name = all_group_array[2]
    }
    
    
    index_control_TMEresults = which(rownames(rt) %in% group$sample[group$group==group_control_name])
    index_case_TMEresults = which(rownames(rt) %in% group$sample[group$group==group_case_name])
    
    
    immune_cell_P.val_list_all = list()
    print('--默认【蓝色→正常组，红色→实验组】--')
    #对每个免疫细胞循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
    total_num = ncol(rt)
    for(i in 1:total_num){
        Pett_ProgressBar(i, total_num)
        
        normalData=rt[index_control_TMEresults, i]
        tumorData=rt[index_case_TMEresults, i]
        #画正常组
        vioplot(normalData, at=3*(i-1),lty=1,add = T,col = 'blue')
        #画肿瘤组
        vioplot(tumorData, at=3*(i-1)+1,lty=1,add = T,col = 'red')
        #非参数检验
        wilcoxTest = wilcox.test(normalData, tumorData)
        #p值 取3位小数
        p=round(wilcoxTest$p.value, 3)
        # print(paste0('--',colnames(rt)[i], '的P值[', p, ']--'))
        immune_cell_P.val_list_all[[colnames(rt)[i]]] = p
        
        #p值下面小横线位置，用每组数据里每种免疫细胞类型的最大占比定y的位置
        mx=max(c(normalData, tumorData))
        #画小短线
        lines(c(x=3*(i-1)+0.2, x=3*(i-1)+0.8), c(mx, mx))
        #标注p值, p<0.001 时标p<0.001， 否则标出p的数值
        text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), 
             cex = 0.8) # p值字体大小，这步在ggplot2中有点难搞……
        #标横坐标不同细胞类型
        # text(seq(1,64,3), -0.05, xpd = NA, labels=colnames(rt), cex = 1, srt = 45, pos=2)
        text(seq(1,64,3), - abscissa_distance, xpd = NA, labels=colnames(rt), cex = 1, srt = 45, pos=2)
    }
    
    
    print('---【图片，请到工作目录找pdf文件】---')
    print(paste0('~~文件名为【', file_name, '】~~'))
    dev.off()
    
    immune_cell_p_sig0.05_list = list()
    # immune_cell_p_sig0.05 = NULL
    immune_cell_name_array = get_list_names(immune_cell_P.val_list_all)
    for(a_cell in immune_cell_name_array){
        p.val = immune_cell_P.val_list_all[[a_cell]]
        if(p.val < 0.05){
            immune_cell_p_sig0.05_list[[a_cell]] = p.val
        }
    }
    return(     list(all=immune_cell_P.val_list_all, p_sig0.05=immune_cell_p_sig0.05_list)       )
}




two_array = function(array1, array2='', action=''){
    
    if( class(array1) != 'list' ){
        array_all = list(array1, array2)
    }else{
        array_all = array1
    }
    
    array_final = 'Pett_15270928527'
    for( array0 in array_all ){
        if( array_final[1] == 'Pett_15270928527' ){
            array_final = array0
            next
        }
        
        if( action=='交' ){
            array_final = intersect(array_final, array0)
        }
        else if( action %in% c('并', '+') ){
            array_final = union(array_final, array0)
        }
        else if( action %in% c('差', '减', '不同', '不', '-')  ){
            print('---[求向量x与向量y中不同的元素(只取x中不同的元素)]---')
            array_final = setdiff(array_final, array0)
        }
        else if(  action %in% c('同', '一样', '相同', '=')  ){
            array_final = setdiff(array_final, array0)
        }else{
            print('---【输入的action有误，或者没有输入action！】---')
        }
    }   ##### 循环结束
    return(array_final)
}




### Plot Figure 4A
immune_cell_P.val_list_list = plot_vioplot_immune( TME.results )

immune_cell_P.val_list_all = immune_cell_P.val_list_list[[1]]
immune_cell_p_sig0.05_list = immune_cell_P.val_list_list[[2]]
(positive_immune_cell_array = get_list_names(immune_cell_p_sig0.05_list))
positive_immune_cell_array = two_array(positive_immune_cell_array, c('B cells naive', 'B cells memory', 'Mast cells resting', 'Monocytes', 'Macrophages M1'), action = '差')
# positive_immune_cell_array = two_array(positive_immune_cell_array, c('B cells naive', 'B cells memory', 'Mast cells resting', 'Monocytes'), action = '差')


TME_results_chosed = TME.results[,   colnames(TME.results) %in% get_list_names(immune_cell_p_sig0.05_list)       ]




Correlation_calculate = function(axis_x, axis_y=NULL, pearson_spearman=c('pearson','spearman')[1]     ){
    # https://www.cnblogs.com/timeisbiggestboss/p/8477138.html
    library('psych') #### 输入数据可以是data.frame，输出两两变量的相关系数R值，显著性水平a值
    
    # library(Hmisc)#加载包
    # 包里的rcorr()函数能够同时给出相关系数以及显著性水平p-value。 rcorr(x, type = c(“pearson”,“spearman”))。
    mydata <- mtcars[, c(1,3,4,5,6,7)]
    if(is.null(axis_y)){
        ### 方法一：用【axis_x】的columns相互进行Correlation计算
        cat('用【axis_x】的columns相互进行Correlation计算')
        Correlation_result_list = corr.test(axis_x, method=pearson_spearman)
    }else{
        ### 方法二：用【axis_x】的columns与【axis_y】的columns进行Correlation计算
        rowname_intersect_array = intersect(  rownames(axis_x), rownames(axis_y)   )
        cat(paste0('两df的sample相交有【', length(rowname_intersect_array), '】个\naxis_x的rownames有【', length(rownames(axis_x)), '】个'    ))
        cat('用【axis_x】的columns与【axis_y】的columns进行Correlation计算')
        Correlation_result_list = corr.test(x=axis_y[ rownames(axis_y) %in% rowname_intersect_array ,], y=axis_x[ rownames(axis_x) %in% rowname_intersect_array ,], method=pearson_spearman)
    }
    
    ### Correlation_result_list$P → P值
    ### Correlation_result_list$r → r相关系数
    return(Correlation_result_list)
}


library(readxl)
Table_3 = read_excel('./data/Supplementary Table 3.xlsx', sheet=1, na=c("NA","#DIV/0!") )
QADGs_array = Table_3$gene_symbol

###上调
Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(QADGs_array[QADGs_array %in% rownames(ADGs)[ADGs$logFC > 0]])  ,]), axis_y=TME_results_chosed[,colnames(TME_results_chosed) %in% positive_immune_cell_array], pearson_spearman=c('pearson','spearman')[1])


plot_Correlation_many_in_one = function(exp_t_bind_TME.results, axis_x_y_list, max_col_num=5 ){
    warning('---【直接运行函数，不能生成本地png】---')
    if(F){
        
        
        max_col_num = 6
        print('exp_t_bind_TME.results = as.data.frame(cbind(  t(exp),    TME.results  ))')
        # Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  c(gene_up, gene_down)  ,]), axis_y=TME_results_chosed, pearson_spearman=c('pearson','spearman')[1])
        # Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(c(drug_ADDEGs_gene_array ))  ,]), axis_y=TME_results_chosed, pearson_spearman=c('pearson','spearman')[1])
        # Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(c(F_ADGs_LC1.2_up_array, F_ADGs_LC1.2_down_array))  ,]), axis_y=TME_results_chosed, pearson_spearman=c('pearson','spearman')[1])
        # Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(lassogene)  ,]), axis_y=TME_results_chosed, pearson_spearman=c('pearson','spearman')[1])
        # 
        # axis_x_y_list = list()
        # positive_immune_cell_array = c( 'B cells memory', 'T cells regulatory (Tregs)', 'NK cells resting', 'Macrophages M2', 'Mast cells resting', 'Mast cells activated', 'Eosinophils',  'Neutrophils')
        # for(   index in 1:length(positive_immune_cell_array)   ){
        #     test1 = Correlation_result_list$r[index  ,]
        #     test2 = Correlation_result_list$p[index  ,]
        #     test2_new = test2[test2 < 0.05]
        #     immune_gene_chosed_array = test1[abs(test1) >= 0.3  &  names(test1) %in% names(test2_new)]
        #     if(length(immune_gene_chosed_array) > 0){
        #         axis_x_y_list[[positive_immune_cell_array[index]]] = names(immune_gene_chosed_array)
        #     }
        # }
        
        library(ggpubr)
        plotlist = list()
        
        if(    class(axis_x_y_list) == 'list'    ){
            axis_x_array = get_list_names(axis_x_y_list)
            total_num = length(axis_x_array)
            
            num_plot = 0
            for(   num in 1:total_num   ){
                Pett_ProgressBar( num,total_num )
                axis_x = axis_x_array[num]
                axis_y_array = axis_x_y_list[[axis_x]]
                
                for(num_y in 1:length(axis_y_array) ){
                    axis_y = axis_y_array[num_y]
                    
                    num_plot = num_plot + 1
                    plotlist[[num_plot]] =  ggscatter(exp_t_bind_TME.results, x = axis_x, y = axis_y,
                                                      color = "blue", size = 1, # Points color and size
                                                      add = "reg.line",  # Add regression line
                                                      add.params = list(color = "red", fill = "gray"), # Customize regression line
                                                      conf.int = TRUE, # Add confidence interval
                                                      cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                                      cor.coeff.args = list(method = "pearson"))
                }
                
                
                
            }
            
            
        }
        
        
        
        library('cowplot')
        file_name = paste0("相关性分析-多个[", get_date(), ' ', get_time(for_write_file_in_windows=T),"].png")
        # pdf(  file_name,  height=8, width=10   )              #保存图片的文件名称
        # png(file=file_name, height=1000, width=1000, bg="transparent")
        # png(file=file_name, height=600, width=1200, bg="white")
        max_col_num = ifelse(  max_col_num >= length(plotlist), length(plotlist),  max_col_num    )
        num_of_row = ceiling(   length(plotlist) / max_col_num   )
        height = num_of_row * 200
        width = ifelse(  length(plotlist) > max_col_num, 5 * 230,  length(plotlist) * 230   )
        
        png(file=file_name, height=height, width=width, bg="white") # white, transparent#
        png_plot_grid = plot_grid(plotlist = plotlist,
                                  # labels = as.character(1:length(plotlist)),
                                  # nrow = 1,
                                  ncol = max_col_num,
                                  align = "h"
        )
        print(png_plot_grid)
        dev.off()
        print('---【图片，请到工作目录找pdf文件】---')
        print(paste0('~~文件名为【', file_name, '】~~'))
        warning('---【直接运行函数，不能生成本地png】---')
        
        
        # return(plotlist)
    }
}

x_y_immune_list = list()
for(   index in 1:length(positive_immune_cell_array)   ){
    Pett_ProgressBar(index, length(positive_immune_cell_array))
    test1 = Correlation_result_list$r[index  ,]
    test2 = Correlation_result_list$p[index  ,]
    test2_new = test2[test2 < 0.05]
    immune_gene_chosed_array = test1[abs(test1) >= 0.3  &  names(test1) %in% names(test2_new)]
    if(length(immune_gene_chosed_array) > 0){
        x_y_immune_list[[positive_immune_cell_array[index]]] = names(immune_gene_chosed_array)
    }
}
x_y_immune_list




plot_Correlation_multi_in_one = function(   Correlation_result_list, axis_x_y_list='', title='', cut_with_P.val_level=0.05, cut_with_R.val_level=0.03, order=c("original","AOE", "FPC", "hclust", "alphabet")[1]  , height=10, width=20    ){
    if(  axis_x_y_list != ''   &  class(axis_x_y_list) %in% c('list')   ){
        Corr_r = Correlation_result_list$r
        Corr_p = Correlation_result_list$p
        if(   length(intersect(rownames(Corr_p), get_list_names(axis_x_y_list)))  > 0   ){
            for(  a_row in rownames(Corr_p)  ){
                cat(a_row)
                cat('\n')
                Corr_p[   a_row   ,   ! colnames(Corr_p) %in%  axis_x_y_list[[a_row]]   ] = 100
            }
            
            Eligible_columns_index1 = ! colSums(Corr_p >= cut_with_P.val_level) == dim(Corr_p)[1]
            Eligible_gene1 = colnames(Corr_p)[Eligible_columns_index1]
            Eligible_columns_index2 = ! colSums(Corr_r <= cut_with_R.val_level) == dim(Corr_r)[1]
            Eligible_gene2 = colnames(Corr_r)[Eligible_columns_index2]
            
            sig_gene_array = intersect(Eligible_gene1, Eligible_gene2)
            Corr_r = Corr_r[,  colnames(Corr_r)  %in%  sig_gene_array  ]
            Corr_p = Corr_p[,  colnames(Corr_p)  %in%  sig_gene_array  ]
            
            ### 重新生成【Correlation_result_list】
            Correlation_result_list = list(r=Corr_r, p=Corr_p)
            
        }else{
            cat(paste0('\n', '--【项目未完成】--'))
            for(  a_col in colnames(Corr_p)  ){
                cat(a_col)
                cat('\n')
                Corr_p[    ! rownames(Corr_p) %in%  axis_x_y_list[[a_col]]  ,     a_col   ] = 100
            }
        }
    }
    
    
    
    library(corrplot)
    print('--【需要先用[Correlation_calculate]函数生成Correlation_result_list】--')
    # order参数设定不同展示顺序，默认order="orginal"以原始顺序展示
    
    file_name = paste0("Correlation[", get_date(), ' ', get_time(for_write_file_in_windows=T),"].pdf")
    pdf(  file_name,  height=10, width=20)              #保存图片的文件名称
    
    
    if(   is.numeric(cut_with_P.val_level)   ){
        print(paste0('--结合P值【', cut_with_P.val_level,'】绘制 相关图--'))
        png_Correlation = corrplot(Correlation_result_list$r, title=title, order=order, p.mat = Correlation_result_list$p, sig.level = cut_with_P.val_level, insig = "blank") # method = "number"
    }else{
        print(paste0('--【不结合】P值【', cut_with_P.val_level,'】绘制 相关图--'))
        png_Correlation = corrplot(Correlation_result_list$r, title=title, order=order, p.mat = Correlation_result_list$p,   insig = "blank")
    }
    
    print('---【图片，请到工作目录找pdf文件】---')
    print(paste0('~~文件名为【', file_name, '】~~'))
    dev.off()
    
    return(png_Correlation)
}

### Plot Figure 4B
png_correlation = plot_Correlation_multi_in_one(   Correlation_result_list, x_y_immune_list, title='',  height=5, width=12,  cut_with_P.val_level=0.05, cut_with_R.val_level=0.03, order=c("original","AOE", "FPC", "hclust", "alphabet")[1]    )



###下调
Correlation_result_down_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(QADGs_array[QADGs_array %in% rownames(ADGs)[ADGs$logFC < 0]])  ,]), axis_y=TME_results_chosed[,colnames(TME_results_chosed) %in% positive_immune_cell_array], pearson_spearman=c('pearson','spearman')[1])


x_y_immune_list = list()
for(   index in 1:length(positive_immune_cell_array)   ){
    Pett_ProgressBar(index, length(positive_immune_cell_array))
    test1 = Correlation_result_list$r[index  ,]
    test2 = Correlation_result_list$p[index  ,]
    test2_new = test2[test2 < 0.05]
    immune_gene_chosed_array = test1[abs(test1) >= 0.3  &  names(test1) %in% names(test2_new)]
    if(length(immune_gene_chosed_array) > 0){
        x_y_immune_list[[positive_immune_cell_array[index]]] = names(immune_gene_chosed_array)
    }
}
x_y_immune_list


### Plot Figure 4C
png_correlation = plot_Correlation_multi_in_one(   Correlation_result_list, x_y_immune_list, title='',  height=5, width=12,  cut_with_P.val_level=0.05, cut_with_R.val_level=0.03, order=c("original","AOE", "FPC", "hclust", "alphabet")[1]    )






CortexRelated_QADGs_array = c('GSR', 'MAT2A', 'ITGA8', 'SLC2A1', 'DLG5', 'MYOD1', 'ID2', 'CAT', 'PANX1', 'DNAJC3', 'KLHL21', 'MLEC', 'NACC2', 'TNRC6A', 'KNOP1', 'FSTL1', 'LIX1', 'ITPRIPL2', 'IRF2BPL', 'RAPGEF5', 'VCAN', 'FOXO1', 'TBC1D2B', 'SASH1', 'PTGR1', 'CACHD1', 'COL11A1', 'DDIT4L', 'GPR37', 'IL4R', 'TMEM132B', 'HSP90B1', 'TFE3', 'NR2E1', 'KANK1', 'ST5', 'PLA2G4A', 'TLE4', 'ZEB2', 'PAX6', 'FAM167A', 'GPER1', 'ZSWIM6', 'LGR4', 'METTL7A', 'DNAJB9', 'ADRBK1', 'PPP5C', 'GPHN', 'ZMYND11', 'FXR1', 'INPP5A', 'EXOC2', 'LAMA5', 'ZFP91', 'SMARCD3', 'RNF41', 'PALD1', 'AP2A1', 'NOTCH3', 'IDH3A', 'EEF2', 'HDGF', 'CXCL16', 'RNF208')
DEGs_chosed_ENTREZID_array3 = gene_convert_symbol_ENTREZID(unique( CortexRelated_QADGs_array  ))$ENTREZID
KEGG_p.val_cutoff = 1
KEGG_q.val_cutoff = 1
organism = ifelse(tolower(Taxonomy) %in% c('homo', 'homo sapiens', 'human', '9606', 9606), "hsa", 'mmu')


ego <- enrichGO(
    # gene = DEGs_chosed_ENTREZID_array[!is.na(DEGs_chosed_ENTREZID_array)],
    # gene = all_gene_ENTREZID_array,
    # gene = all_input_gene_ENTREZID_array ,
    gene = DEGs_chosed_ENTREZID_array3 ,
    # universe      = DEG_symbol_ENTREZID$ENTREZID,
    OrgDb         = OrgDb,
    ont           = toupper(BP_MF_CC_ALL), ### One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    # ont           = toupper('ALL'), ### One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    ### BP（生物学过程）、MF（分子功能）和CC（细胞组成）
    ### 
    ### 分子功能(Molecular Function,MF )
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_p.val_cutoff,
    qvalueCutoff  = GO_q.val_cutoff,
    readable      = TRUE)

png_GO_barplot = barplot(ego, title = "Biological process", showCategory=20)

library('cowplot')
png(file='Figure 5C.png', height=600, width=400, bg="white")
print(png_GO_barplot)
dev.off()





Undetermined_QADGs_array = c('CDK1', 'LBP', 'CSK', 'TGFB1', 'HGF', 'BLNK', 'LRP1', 'CDK7', 'ESPL1', 'MCM4', 'KLKB1', 'ADH4', 'NUF2', 'CDCA5', 'NINL', 'TPX2', 'MLC1', 'ELMO1', 'SEPT1', 'AMER1', 'EPHB3', 'CAV1', 'CASP7', 'XIAP', 'FOXO6', 'FGR', 'S1PR3', 'GPX3', 'STAT4', 'PRC1', 'CEBPD', 'IL1R1', 'CXCR4', 'LAMB2', 'TRPV4', 'STMN1', 'BDNF', 'CDK5', 'ACTN2', 'VEGFA', 'UHRF1', 'CCNA1', 'AREG', 'SRA1', 'UBE2H', 'UBE2T', 'USF1', 'IRF7', 'KLF6', 'LSM4')
DEGs_chosed_ENTREZID_array4 = gene_convert_symbol_ENTREZID(unique( Undetermined_QADGs_array  ))$ENTREZID
KEGG_p.val_cutoff = 1
KEGG_q.val_cutoff = 1
organism = ifelse(tolower(Taxonomy) %in% c('homo', 'homo sapiens', 'human', '9606', 9606), "hsa", 'mmu')

print('Make sure the Internet is connected!')
kegg_data2 <- enrichKEGG(
    # gene = DEGs_chosed_ENTREZID_array[!is.na(DEGs_chosed_ENTREZID_array)], ### gene = gene_up$ENTREZID,
    # gene = all_gene_ENTREZID_array,
    gene = DEGs_chosed_ENTREZID_array4,
    # universe      = DEG_symbol_ENTREZID$ENTREZID,
    # gene = symbol_ENTREZID$ENTREZID, ### gene = gene_up$ENTREZID,
    organism = organism, # homo→hsa；mice→mmu；
    keyType = "kegg", # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
    # universe = gene_all,
    pvalueCutoff = KEGG_p.val_cutoff,
    qvalueCutoff = KEGG_q.val_cutoff) #P值或者Q值阈值可以根据富集的结果进行主观修改

png_KEGG_dotplot2 = enrichplot::dotplot( kegg_data2, title = "KEGG pathways", showCategory = 20, )#气泡图
png(file='Figure 6B.png', height=600, width=400, bg="white")
print(png_KEGG_dotplot2)
dev.off()












##### MMSE #####
pheno_frontal = pheno
# pheno_frontal$MMSE = as.numeric(pheno_frontal$`mmse:ch1`)
pheno_frontal$MMSE = as.numeric(pheno_frontal$`MMSE`)
pheno_frontal = pheno_frontal[ !is.na(pheno_frontal$MMSE) ,]
pheno_frontal = pheno_frontal[  pheno_frontal$geo_accession %in% group$sample  ,]


pheno_frontal_10_26 = subset(pheno_frontal,  MMSE >= 10   &    MMSE <= 26  )

pheno_MMSE_score_chosed = data.frame(MMSE=as.numeric(pheno_frontal_10_26$MMSE))
rownames(pheno_MMSE_score_chosed) = pheno_frontal_10_26$geo_accession


immuneRelated_QADGs = c(  'CDC37', 'RRAGB', 'WNT2', 'AIP', 'THRB', 'WASL', 'NEK9', 'PLCG1', 'RAP1B', 'RORA', 'ALDH2', 'CD81', 'STAT3', 'ERBB3', 'EPS8', 'MMP16', 'RAC3', 'SMAD3', 'TF', 'FYN', 'FXR1', 'MYOF', 'SGK1', 'MYO6', 'PTGDS', 'ABL1', 'RRAGC', 'CDK8', 'DEK', 'NCOR2', 'XPO1', 'FOXO3', 'KIFC1', 'GGA3', 'ANXA5', 'MED1', 'ARNT2', 'ATF4', 'HMGB1', 'TEX10', 'CWC27', 'RNPC3', 'CRYAB', 'RAF1', 'RHOB', 'SMC2', 'HSPB3', 'CRKL', 'DLG3', 'RAD1', 'NEK2', 'GRIA4', 'NGF', 'HDAC9', 'LRRC7', 'TRIP6', 'STS', 'HSPB1', 'CCNB1', 'VWF', 'RBL2', 'GPC1', 'FGF2', 'ANLN', 'LRIG1', 'NTRK2', 'TUB', 'GRM5', 'SOX2', 'NOS2', 'EPHA3', 'RTN4', 'PLCB1', 'MEF2C', 'SDC3', 'PLS1', 'RHOU', 'NEK7', 'EFNA1', 'BAG3', 'SDC4', 'HDAC1', 'LPAR1', 'HIF3A', 'EDN1', 'FOXO1', 'RFC3', 'AKAP5', 'FHL2', 'DCC', 'SYNJ1', 'FGFR2', 'DAPP1', 'ASAP1', 'GRIA3', 'WASF2', 'EGR2', 'KAT2B', 'BAG4', 'IDH3A', 'SRRM2', 'WDR5', 'MAEL', 'PREP', 'GREB1', 'AHSA2', 'GRM2', 'SMYD3', 'AAAS', 'BRD1', 'MYOD1', 'LAMA5', 'EPAS1', 'NSUN2', 'GRIA2', 'UBB', 'S1PR1', 'ARF4', 'MCM7', 'TCHP', 'CDKN3', 'FGFR3', 'LBR', 'GIT2', 'ACTB', 'COX2', 'UTRN', 'PTPRK', 'IRS2', 'DVL1', 'SAT1', 'LAMB3', 'ACTN1', 'NRG4', 'BCL6', 'PRR5', 'RAB1A', 'LPAR3', 'PTPRJ', 'PDGFC', 'HIF1A', 'DAG1', 'PCGF5', 'CEBPB', 'BRD8', 'SNW1', 'NOP58', 'UBE2S', 'FLCN', 'ATP6', 'SET', 'UBE2N', 'EIF3L', 'CALM1', 'ATP5E'  )



for(input_gene_str in c('Immune_Related_QADGs', 'Cortex_Related_QADGs', 'Undetermined_QADGs')){


if(input_gene_str == 'Cortex_Related_QADGs'){
    # input_gene_array = gene06_CortexRelated_QADGs
    input_gene_array = CortexRelated_QADGs_array
    # Pett_read_write_table('CortexRelated_QADGs.xlsx', CortexRelated_QADGs)
}else if(input_gene_str == 'Immune_Related_QADGs'){
    # input_gene_array = gene05_Immune_Related_QADGs
    #       ↓↓↓↓↓↓QADGs_array还是lassogene
    input_gene_array = immuneRelated_QADGs
    # Pett_read_write_table('immuneRelated_QADGs.xlsx', immuneRelated_QADGs)
    
}else if(input_gene_str == 'Undetermined_QADGs'){
    input_gene_array = Undetermined_QADGs_array
    # Pett_read_write_table('Undetermined_QADGs.xlsx', Undetermined_QADGs)
}



Correlation_result_list = Correlation_calculate(axis_x=t(exp[  rownames(exp) %in%  unique(input_gene_array)  ,   colnames(exp) %in% rownames(pheno_MMSE_score_chosed)  ]), axis_y=pheno_MMSE_score_chosed, pearson_spearman=c('pearson','spearman')[1])






positive_correlation_MMSE_gene_array = NULL
MMSE_array = c('MMSE')
x_y_MMSE_list = list()
for(   index in 1:length(MMSE_array)   ){
    test1 = Correlation_result_list$r[index  ,]
    test2 = Correlation_result_list$p[index  ,]
    test2_new = test2[test2 < 0.05]
    immune_gene_chosed_array = test1[abs(test1) >= 0.3  &  names(test1) %in% names(test2_new)]
    if(length(immune_gene_chosed_array) > 0){
        x_y_MMSE_list[[MMSE_array[index]]] = names(immune_gene_chosed_array)
        positive_correlation_MMSE_gene_array = c(positive_correlation_MMSE_gene_array, names(immune_gene_chosed_array))
    }
}
x_y_MMSE_list
positive_correlation_MMSE_gene_array
cat(paste0(positive_correlation_MMSE_gene_array, collapse = '\n'))




exp_t_bind_TME.results = as.data.frame(cbind(  t(exp[, colnames(exp)    %in%  rownames(pheno_MMSE_score_chosed)   ]),    pheno_MMSE_score_chosed  ))

axis_x_y_list = x_y_MMSE_list

max_col_num = 6

library(ggpubr)
plotlist = list()

if(    class(axis_x_y_list) == 'list'    ){
    axis_x_array = get_list_names(axis_x_y_list)
    total_num = length(axis_x_array)
    
    num_plot = 0
    for(   num in 1:total_num   ){
        Pett_ProgressBar( num,total_num )
        axis_x = axis_x_array[num]
        axis_y_array = axis_x_y_list[[axis_x]]
        
        for(num_y in 1:length(axis_y_array) ){
            axis_y = axis_y_array[num_y]
            
            num_plot = num_plot + 1
            plotlist[[num_plot]] =  ggscatter(exp_t_bind_TME.results, x = axis_x, y = axis_y,
                                              color = "blue", size = 1, # Points color and size
                                              add = "reg.line",  # Add regression line
                                              add.params = list(color = "red", fill = "gray"), # Customize regression line
                                              conf.int = TRUE, # Add confidence interval
                                              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                              cor.coeff.args = list(method = "pearson"))
        }
        
        
        
    }
    
    
}



library('cowplot')
file_name = paste0( './data/Figure 7 - ', input_gene_str, '.png' )
# pdf(  file_name,  height=8, width=10   )              #保存图片的文件名称
# png(file=file_name, height=1000, width=1000, bg="transparent")
# png(file=file_name, height=600, width=1200, bg="white")
max_col_num = ifelse(  max_col_num >= length(plotlist), length(plotlist),  max_col_num    )
num_of_row = ceiling(   length(plotlist) / max_col_num   )
height = num_of_row * 200
width = ifelse(  length(plotlist) > max_col_num, 5 * 230,  length(plotlist) * 230   )

png(file=file_name, height=height, width=width, bg="white") # white, transparent#
png_plot_grid = plot_grid(plotlist = plotlist,
                          # labels = as.character(1:length(plotlist)),
                          # nrow = 1,
                          ncol = max_col_num,
                          align = "h"
)
png_plot_grid
dev.off()
print('---【图片，请到工作目录找pdf文件】---')
print(paste0('~~文件名为【', file_name, '】~~'))





}










##### miRNA #####
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(stringr))
exprSet <- read.table( './data/GSE48552_normalized_counts.txt.gz', header = T, row.names = 1 )

eSet_GSE48552 = getGEO( 'GSE48552', destdir = './data',  # downloaded in [#02Annotation] folder
                       AnnotGPL = F, getGPL = F)  # remove platform info


save( eSet_GSE48552, file = './data/eSet_eSet_GSE48552.rdata')



load('./data/eSet_eSet_GSE48552.rdata')
pheno_GSE48552 = eSet_GSE48552[[1]]@phenoData@data


group_df = pheno_GSE48552[, c('geo_accession', 'source_name_ch1')]
group_df$batch_num = 1
colnames(group_df)[c(1,2)] = c('sample', 'group')
group_df$group = gsub(' ', '_', gsub('_prefrontal cortex', '', group_df$group))



suppressMessages(library('DESeq2'))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library("amap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("BiocParallel"))
# suppressMessages(library("YSX"))
suppressMessages(library("ImageGP"))
suppressMessages(library('sva'))
suppressMessages(library('ggfortify'))
suppressMessages(library('patchwork'))
suppressMessages(library('ggbeeswarm'))



Pett_deseq2normalizedExpr <- function(dds, output_prefix='ehbio', rlog=T, vst=F, savemat=F){
    #标准化后的结果按整体差异大小排序，同时输出对数转换的结果。
    # message( "rlog() may take a few minutes with 30 or more samples，vst() is a much faster transformation" )
    message("rlog() may take a long time with 50 or more samples，vst() is a much faster transformation")
    cat('---【数据集小于30 -> rlog，大数据集大于30 -> vst】---')
    message('此函数做了log2(exp(1))——https://www.jianshu.com/p/cd7aa2d6b77b')
    
    # Get normalized counts
    normalized_counts <- DESeq2::counts(dds, normalized=TRUE) ## “FALSE”意思指各处理之间的差异不增加试验的方差平均数
    
    # 标准化的结果按整体差异大小排序
    normalized_counts_mad <- apply(normalized_counts, 1, mad)
    normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
    
    # 常规R输出忽略左上角（输出文件第一列也就是基因列的列名字）
    # 输出结果为完整矩阵，保留左上角的id。
    ##### 【# 标准化后的数据输出】 #####
    normalized_counts_output = data.frame(id=rownames(normalized_counts), normalized_counts)
    
    if(savemat){
        print("Output normalized counts")
        write.table(normalized_counts_output, file=paste0(output_prefix,".DESeq2.normalized.xls"),
                    quote=F, sep="\t", row.names=F, col.names=T)
    }
    
    normexpr <- list(normalized=normalized_counts, normalizedSave=normalized_counts_output)
    
    if(rlog) {
        rld <- DESeq2::rlog(dds, blind=FALSE)
        rlogMat <- assay(rld)
        rlogMat_mad <- apply(rlogMat, 1, mad)
        rlogMat <- rlogMat[order(rlogMat_mad, decreasing=T), ]
        
        rlogMat_output = data.frame(id=rownames(rlogMat), rlogMat)
        
        if(savemat) {
            print("Output rlog transformed normalized counts")
            write.table(rlogMat_output, file=paste0(output_prefix,".DESeq2.normalized.rlog.xls"),
                        quote=F, sep="\t", row.names=F, col.names=T)
        }
        normexpr$rlog <- rlogMat
        normexpr$rlogSave <- rlogMat_output
        
    }
    
    
    if( vst ) {
        rld <- DESEq2::vst(dds, blind=FALSE)
        vstMat <- assay(rld)
        vstMat_mad <- apply(vstMat, 1, mad)
        vstMat <- vstMat[order(vstMat_mad, decreasing=T), ]
        
        vstMat_output = data.frame(id=rownames(vstMat), vstMat)
        
        if(savemat){
            print("Output vst transformed normalized counts")
            write.table(vstMat_output, file=paste0(output_prefix,".DESeq2.normalized.vst.xls"),
                        quote=F, sep="\t", row.names=F, col.names=T)
        }
        
        normexpr$vst <- vstMat
        normexpr$vstSave <- vstMat_output
    }
    
    
    return(normexpr)
}






Pett_Base_Matrix2colCorrelation <-
    function(mat,
             method = "pearson",
             digits = 4,
             cor_file = NULL) {
        pearson_cor <-
            round(as.matrix(cor(mat, method = method)), digits = digits)
        hc <- amap::hcluster(t(mat), method = method)
        pearson_cor <- pearson_cor[hc$order, hc$order]
        if (!is.null(file)  &  FALSE ) { ##### 不输出文件
            pearson_cor_output = data.frame(id = rownames(pearson_cor), pearson_cor)
            write.table(
                pearson_cor_output,
                file = cor_file,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = T
            )
        }
        return(list(pearson_cor = pearson_cor, hc = hc))
    }








Plot_clusterSampleHeatmap2 <- function (mat, method="pearson", digits=4,
                                        cor_file=NULL, saveplot=NULL, ...){
    print("Performing sample clustering")
    
    hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(100)
    
    mat <- dataFilter(mat, ...)
    
    pearson_cor_hc <- Pett_Base_Matrix2colCorrelation(mat, method, digits, cor_file)
    pearson_cor <-  pearson_cor_hc$pearson_cor
    hc <-  pearson_cor_hc$hc
    
    if(!is.null(saveplot)) {
        base_plot_save(saveplot, ...)
    }
    
    gplots::heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
                      col=hmcol, margins=c(11,11), key=T,
                      main="Correlation plot")
    if(!is.null(saveplot)) {
        dev.off()
    }
}






Plot_clusterSampleUpperTriPlot <- function (mat, method="pearson", digits=4,
                                            cor_file=NULL, saveplot=NULL,
                                            width=13.5, height=15, ...){
    print("Performing sample clustering")
    
    pearson_cor_hc <- Pett_Base_Matrix2colCorrelation(mat, method, digits, cor_file)
    pearson_cor <-  pearson_cor_hc$pearson_cor
    
    upper_tri <- get_upper_tri(pearson_cor)
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
    
    col = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(100)
    
    # Create a ggheatmap
    p <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradientn(colours=col, name=paste(method,"correlation")) +
        theme_classic() +
        coord_fixed() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.justification = c(1, 0),
            legend.position = "top",
            legend.direction = "horizontal")+
        guides(fill = guide_colorbar(barwidth = 9, barheight = 1,
                                     title.position = "left"))
    
    if(!is.null(saveplot)){
        ggsave(plot=p, filename=saveplot, units=c("cm"),...)
    }
    return(p)
}



# count_matrix_file = 'ehbio_trans.Count_matrix.txt'
# sampleFile = 'sampleFile.txt'



Pett_readscount2deseq <- function(count_matrix_file, sampleFile, design, covariate=NULL,
                                  filter=F, rundeseq=T ) {
    if( ! class(count_matrix_file) %in% 'data.frame'  ){
        data <- read.table(count_matrix_file, header=T, row.names=1, com='', quote='',
                           check.names=F, sep="\t")
    }else{
        data = count_matrix_file
        tryCatch({       rm(count_matrix_file)        }, error=function(cond){        1111111       })
    }
    
    
    
    if( ! class(sampleFile) %in% 'data.frame'  ){
        sample <- read.table(sampleFile, header=T, row.names=1, com='',
                             quote='', check.names=F, sep="\t")
        sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
    }else{
        sample = sampleFile
        tryCatch({       rm(sampleFile)        }, error=function(cond){        1111111       })
    }
    
    
    
    if(  !is.null(covariate)  ){
        covariate <- paste(covariate, collapse="+")
        formula <- as.formula(paste("~", covariate,"+", design))
    } else {
        formula <- as.formula(paste("~", design))
    }
    
    
    
    dds <- DESeqDataSetFromMatrix(countData = data, colData = sample, design=formula)
    
    print(paste("Read in", nrow(dds),"genes"))
    
    if(  is.null(filter)  ){
        filter = nrow(sample)/2
    }
    
    if(  is.numeric(filter)  ){
        keep <- rowSums(DESeq2::counts(dds)) > filter
        dds <- dds[keep,]
        print(paste('【', nrow(dds),"】genes remained after filtering of genes with all counts less than", nrow(sample)/2, "in all samples."))
    }
    
    print(paste('【', nrow(dds),"】genes; 【", nrow(sample), "】samples."))
    
    if(  rundeseq  ){
        print("Perform DESeq on given datasets.")
        dds <- DESeq(dds)
    }
    return(dds)
}






DESeqDataSetFromMatrix <- function( countData, colData, design, tidy=FALSE, ignoreRank=FALSE, ... )
{
    
    if (tidy) {
        stopifnot(ncol(countData) > 1)
        rownms <- as.character(countData[,1])
        countData <- countData[,-1,drop=FALSE]
        rownames(countData) <- rownms
    }
    
    # check that these agree in number
    stopifnot(ncol(countData) == nrow(colData))
    
    if (is(countData, "data.frame")) {
        if (any(sapply(countData, is, "factor"))) {
            warning("\n\n  'countData' is a data.frame with one or more factor columns.
  Be aware that converting directly to numeric can lead to errors.
  Provide matrices of integers to avoid this error.")
        }
    }
    
    # we expect a matrix of counts, which are non-negative integers
    countData <- as.matrix( countData )
    
    if (is(colData,"data.frame"))
        colData <- as(colData, "DataFrame")
    
    # check if the rownames of colData are simply in different order
    # than the colnames of the countData, if so throw an error
    # as the user probably should investigate what's wrong
    if (!is.null(rownames(colData)) & !is.null(colnames(countData))) {
        if (all(sort(rownames(colData)) == sort(colnames(countData)))) {
            if (!all(rownames(colData) == colnames(countData))) {
                stop(paste("rownames of the colData:
  ",paste(rownames(colData),collapse=","),"
  are not in the same order as the colnames of the countData:
  ",paste(colnames(countData),collapse=",")))
            }
        }
    }
    if (is.null(rownames(colData)) & !is.null(colnames(countData))) {
        rownames(colData) <- colnames(countData)
    }
    
    se <- SummarizedExperiment(assays = SimpleList(counts=countData), colData = colData, ...)
    object <- DESeqDataSet(se, design = design, ignoreRank)
    
    return(object)
}













title = NULL
scale = TRUE
manual_color_vector = NULL
log_transform = NULL
facet = NULL
size_variable = NULL
color_variable_order = NULL
top_n = 0
shape_variable_order = NULL
dimensions = 2
alpha = 1
label_points = FALSE
log_add = 0
label_font_size = NULL
minimum_mad = 0
debug = FALSE
filename = NULL
legend.position = "right"
coord_fixed_ratio=1


Pett_sp_pca <- function(data,
                        group_data = NULL,
                        title = NULL,
                        scale = TRUE,
                        color_variable = NULL,
                        manual_color_vector = NULL,
                        log_transform = NULL,
                        facet = NULL,
                        size_variable = NULL,
                        shape_variable = NULL,
                        color_variable_order = NULL,
                        top_n = 0,
                        shape_variable_order = NULL,
                        dimensions = 2,
                        alpha = 1,
                        label_points = FALSE,
                        log_add = 0,
                        label_font_size = NULL,
                        minimum_mad = 0,
                        debug = FALSE,
                        filename = NULL,
                        legend.position = "right",
                        coord_fixed_ratio=1,
                        ...) {
    library(ggfortify)
    library(ggplot2)
    
    if (debug) {
        argg <- c(as.list(environment()), list(...))
        print(argg)
    }
    
    if( class(data) == "character" ) {
        data <- sp_readTable(data, row.names = NULL)
        rownames_data <- make.unique(as.vector(data[, 1]))
        data <- data[, -1, drop = F]
        rownames(data) <- rownames_data
    }
    
    data <- data[var(data) != 0, ]
    
    if(  minimum_mad + top_n != 0  ){
        data$mad <- apply(data, 1, mad)
        if(minimum_mad>0){
            data <- data[data$mad > minimum_mad , ]
        }
        data_row_num <- dim(data)[1]
        if (top_n != 0 & top_n < data_row_num) {
            data <-
                data[order(data$mad, decreasing = T),,drop=F]
            data <- data[1:top_n, ]
        }
        data <- data[,-dim(data)[2],drop=F]
    }
    
    data <- as.data.frame(t(data))
    
    #print(data[1:3,1:5])
    
    if(  !sp.is.null(log_transform)  ) {
        # print(y_add)
        # Give the minimum non-zero value to add to avoid log2(0)
        if (log_add == 0) {
            log_add = sp_determine_log_add(data)
        }
        
        data <- data + log_add
        if (log_transform == "log2") {
            data <- log2(data)
        }else if (log_transform == "log10") {
            data <- log10(data)
        }
    }
    
    sampleL = rownames(data)
    
    
    if(  sp.is.null(group_data)  ) {
        data_t_label <- data
        data_t_label$group = sampleL
        data_t_label$Row.names = sampleL
    } else {
        if(  class(group_data) == "character"  ) {
            group_data <- sp_readTable(group_data, row.names = NULL)
            rownames(group_data) <- group_data[, 1]
        }
        #print(colnames(group_data))
        
        data_t_label <- merge(data, group_data, by = 0, all.x = T)
        rownames(data_t_label) <- data_t_label$Row.names
        data_t_label <-
            data_t_label[match(sampleL, data_t_label$Row.names), ]
        
        #print(data_t_label[1:4,1:5])
    }
    
    data_colnames <- colnames(data_t_label)
    
    #print(data_colnames)
    
    if(  !sp.is.null(shape_variable)  ) {
        if(! (shape_variable %in% data_colnames )){
            stop(paste(shape_variable,'must be column names of data!'))
        }
        if (!sp.is.null(shape_variable_order)){
            data_t_label = sp_set_factor_order(data_t_label, shape_variable, shape_variable_order)
        }
        shape_level <- length(unique(data_t_label[[shape_variable]]))
        shapes = (1:shape_level) %% 30
    }
    
    if (!sp.is.null(color_variable)) {
        if(! (color_variable %in% data_colnames )){
            stop(paste(color_variable,'must be column names of data!'))
        }
        if (!sp.is.null(color_variable_order)){
            data_t_label = sp_set_factor_order(data_t_label, color_variable, color_variable_order)
        }
    }
    
    
    pca <- prcomp(data, scale = scale)
    
    rotation = pca$rotation
    x = pca$x
    percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
    
    percentVar2 <- as.data.frame(percentVar)
    rownames(percentVar2) <- colnames(x)
    
    if(  !sp.is.null(filename)  ){
        sp_writeTable(rotation,
                      file = paste0(filename, ".weights.xls"),
                      keep_rownames = T)
        sp_writeTable(x,
                      file = paste0(filename, ".pcs.xls"),
                      keep_rownames = T)
        sp_writeTable(percentVar2,
                      file = paste0(filename, ".pc_variance.xls"),
                      keep_rownames = T)
    }
    
    
    if(  dimensions == 2  ) {
        p = autoplot(pca, data = data_t_label, alpha = alpha, scale=0) + ggtitle(title)
        
        if (!sp.is.null(size_variable)) {
            if(! (size_variable %in% data_colnames )){
                stop(paste(size_variable,'must be column names of data!'))
            }
            size_variable_en = sym(size_variable)
            p <- p + aes(size = !!size_variable_en)
        }
        if (!sp.is.null(color_variable)) {
            color_en = sym(color_variable)
            p <- p + aes(colour = !!color_en)
            p <- sp_manual_color_ggplot2(p,
                                         data,
                                         color_variable,
                                         manual_color_vector)
        }
        
        if (!sp.is.null(shape_variable)) {
            shape_en = sym(shape_variable)
            p <- p + aes(shape = !!shape_en)
            if (shape_level > 6) {
                p <- p + scale_shape_manual(values = shapes)
            }
        }
        
        if (label_points) {
            if (!sp.is.null(label_font_size)) {
                p <-
                    p + geom_text_repel(aes(label = Row.names),
                                        show.legend = F,
                                        size = label_font_size)
            } else {
                p <- p + geom_text_repel(aes(label = Row.names), show.legend = F)
            }
        }
        
        x_label = paste0("PC1 (", round(percentVar[1] * 100,1), "% variance)")
        y_label = paste0("PC2 (", round(percentVar[2] * 100,1), "% variance)")
        
        if(coord_fixed_ratio>0){
            p <- p + coord_fixed(coord_fixed_ratio)
        }
        
        
        p <- sp_ggplot_layout(
            p,
            filename = filename,
            legend.position = legend.position,
            x_label = x_label,
            y_label = y_label,
            title = title,
            ...
        )
        p
        
    } else {
        library(scatterplot3d)
        if (color_variable != "c_t_c_t0304") {
            group = data_t_label[[color_variable]]
            colorA <- rainbow(length(unique(group)))
            
            colors <- colorA[as.factor(group)]
            
            colorl <- colorA[as.factor(unique(group))]
        }
        
        if (shape_variable != "c_t_c_t0304") {
            group <- data_t_label[[shape_variable]]
            pch_l <- as.numeric(as.factor(unique(group)))
            pch <- pch_l[as.factor(group)]
        }
        
        pc <- as.data.frame(pca$x)
        
        saveplot = paste0(filenames, mid, ".pdf")
        
        if (!sp.is.null(saveplot)) {
            base_plot_save(saveplot, ...)
        }
        
        # pdf(paste0(filename,mid,"sds.pdf"))
        scatterplot3d(
            x = pc$PC1,
            y = pc$PC2,
            z = pc$PC3,
            pch = pch,
            color = colors,
            xlab = paste0("PC1 (", round(percentVar[1] * 100), "% variance)"),
            ylab = paste0("PC2 (", round(percentVar[2] * 100), "% variance)"),
            zlab = paste0("PC3 (", round(percentVar[3] * 100), "% variance)")
        )
        
        legend(
            -3,
            8,
            legend = levels(as.factor(color)),
            col = colorl,
            pch = pch_l,
            xpd = T,
            horiz = F,
            ncol = 6
        )
        
        if (!sp.is.null(saveplot)) {
            dev.off()
        }
    }
    
}








Pett_clusterSampleHeatmap2 <- function (mat, method="pearson", digits=4,
                                        cor_file=NULL, saveplot=NULL, ...){
    print("Performing sample clustering")
    
    hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(100)
    
    mat <- dataFilter(mat, ...)
    
    pearson_cor_hc <- Pett_Base_Matrix2colCorrelation(mat, method, digits, cor_file)
    pearson_cor <-  pearson_cor_hc$pearson_cor
    hc <-  pearson_cor_hc$hc
    
    if(!is.null(saveplot)) {
        base_plot_save(saveplot, ...)
    }
    
    gplots::heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
                      col=hmcol, margins=c(11,11), key=T,
                      main="Correlation plot")
    if(!is.null(saveplot)) {
        dev.off()
    }
}






Pett_normalizedExpr2DistribBoxplot <- function(normexpr, saveplot=NULL, ...) {
    if(  'rlog' %in% names(normexpr)  ){
        p <- widedataframe2boxplot(as.data.frame(normexpr$rlog), saveplot=saveplot, ylab="rLog transformed expression value")
    }else if (  'vst' %in% names(normexpr)  ){
        p <- widedataframe2boxplot(as.data.frame(normexpr$vst), saveplot=saveplot, ylab="VST transformed expression value")
    }else if(  'DESeqDataSet' %in% class(normexpr)  ){
        p <- widedataframe2boxplot(as.data.frame(assay(normexpr)), saveplot=saveplot, ylab="VST transformed expression value")
    }
    return(p)
}





comparePairFile=NULL
design="conditions"
padj=0.05
log2FC=1
# "ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
dropCol=c("lfcSE", "stat")
output_prefix="ehbio"


Pett_multipleGroupDEgenes <- function(
    dds,
    comparePairFile=NULL,
    design="conditions",
    padj=0.05,
    log2FC=1,
    # "ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    dropCol=c("lfcSE", "stat"),
    output_prefix="ehbio",
    ...
){
    
    if (file.exists(paste0(output_prefix,".DESeq2.all.DE"))) {
        file.remove(paste0(output_prefix,".DESeq2.all.DE"))
    }
    
    if(  !is.null(comparePairFile)  ){
        compare_data <- read.table(comparePairFile, sep="\t",
                                   check.names=F, quote='', com='')
        colnames(compare_data) <- c("sampA", "sampB")
    } else {
        sampleGroup <- as.data.frame(dds@colData)
        compare_data <- as.vector(unique(sampleGroup[[design]]))
        #compare_data <- letters[1:3]
        len_compare_data <- length(compare_data)
        compareL <- list()
        
        count = 1
        for(i in 1:(len_compare_data-1)) {
            for(j in (i+1):len_compare_data) {
                tmp_compare <- list(sampA=compare_data[i],sampB=compare_data[j])
                compareL[[paste(i,j)]] <- tmp_compare
                count = count + 1
            }
        }
        compare_data <- as.data.frame(do.call(rbind, compareL))
    }
    
    unused <- by(compare_data, 1:nrow(compare_data), function (x)
        twoGroupDEgenes(dds, groupA=unlist(x[1,1]), groupB=unlist(x[1,2]), design=design, padj=padj,
                        log2FC=log2FC, dropCol=dropCol,
                        output_prefix=output_prefix, ...))
    
    #twoGroupDEgenes(dds, tmp_compare, design=design, padj=padj, log2FC=log2FC,
    #                dropCol=dropCol, output_prefix=output_prefix, ...)
}






design="conditions"
padj=0.05
log2FC=1
# "ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
dropCol=c("lfcSE", "stat")
output_prefix="ehbio"



Pett_twoGroupDEgenes <- function
(
    dds,
    groupA,
    groupB,
    design="conditions",
    padj=0.05,
    log2FC=1,
    # "ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    dropCol=c("lfcSE", "stat"),
    output_prefix="ehbio",
    ...
){
    #print(sampleV)
    #groupA <- as.vector(sampleV$sampA)
    #groupB <- as.vector(sampleV$sampB)
    #groupA = 'trt'
    #groupB = 'untrt'
    #design = "conditions"
    #padj = 0.05
    #log2FC = 1
    # "ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    #dropCol = c("lfcSE", "stat")
    #output_prefix = "ehbio"
    
    # print(groupA)
    # if(is.list(groupA)){
    #   groupA <- unlist(groupA)
    # }
    # print(groupA)
    # print(groupB)
    # if(is.list(groupB)){
    #   groupB <- unlist(groupB)
    # }
    # print(groupB)
    print(paste("DE genes between", groupA, groupB, sep=" "))
    contrastV <- c(design, groupA, groupB)
    print(contrastV)
    res <- DESeq2::results(dds,  contrast=contrastV)
    
    normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
    
    
    baseA <- normalized_counts[, colData(dds)[[design]] == groupA]
    if (is.vector(baseA)){
        baseMeanA <- as.data.frame(baseA)
    } else {
        baseMeanA <- as.data.frame(rowMeans(baseA))
    }
    baseMeanA <- round(baseMeanA, 3)
    colnames(baseMeanA) <- groupA
    
    
    baseB <- normalized_counts[, colData(dds)[[design]] == groupB]
    if (is.vector(baseB)){
        baseMeanB <- as.data.frame(baseB)
    } else {
        baseMeanB <- as.data.frame(rowMeans(baseB))
    }
    baseMeanB <- round(baseMeanB, 3)
    colnames(baseMeanB) <- groupB
    
    
    res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
    res <- data.frame(ID=rownames(res), res)
    res$baseMean <- round(rowMeans(cbind(baseA, baseB)),3)
    res$padj[is.na(res$padj)] <- 1
    res$pvalue[is.na(res$pvalue)] <- 1
    res$log2FoldChange <- round(res$log2FoldChange,3)
    res$padj <- as.numeric(formatC(res$padj))
    res$pvalue <- as.numeric(formatC(res$pvalue))
    
    res <- res[order(res$padj),]
    
    comp314 <- paste(groupA, "_vs_", groupB, sep=".")
    
    file_base <- paste(output_prefix, "DESeq2", comp314, sep=".")
    file_base1 <- paste(file_base, "results.xls", sep=".")
    
    #######################################################
    res_output <- res[, !(names(res) %in% dropCol), drop=F]
    #######################################################
    write.table(res_output, file=file_base1, sep="\t", quote=F, row.names=F)
    
    res_de <- res[res$padj<padj, !(names(res) %in% dropCol), drop=F]
    res_de_up <- subset(res_de, res_de$log2FoldChange>=log2FC)
    
    file <- paste(output_prefix, "DESeq2",groupA, "_higherThan_", groupB,
                  'xls', sep=".")
    write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
    
    res_de_up_id <- subset(res_de_up, select=c("ID"))
    file <- paste(output_prefix, "DESeq2",groupA, "_higherThan_", groupB,
                  'id.xls', sep=".")
    write.table(as.data.frame(res_de_up_id), file=file, sep="\t",
                quote=F, row.names=F, col.names=F)
    
    if(dim(res_de_up_id)[1]>0) {
        res_de_up_id_l <- cbind(res_de_up_id, paste(groupA, "_higherThan_",groupB, sep="."))
        write.table(as.data.frame(res_de_up_id_l),
                    file=paste0(output_prefix,".DESeq2.all.DE"),
                    sep="\t",quote=F, row.names=F, col.names=F, append=T)
    }
    
    res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*log2FC)
    file <- paste(output_prefix, "DESeq2",groupA, "_lowerThan_", groupB,
                  'xls', sep=".")
    write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
    
    res_de_dw_id <- subset(res_de_dw, select=c("ID"))
    file <- paste(output_prefix, "DESeq2",groupA, "_lowerThan_", groupB,
                  'id.xls', sep=".")
    write.table(as.data.frame(res_de_dw_id), file=file, sep="\t",
                quote=F, row.names=F, col.names=F)
    
    if(dim(res_de_dw_id)[1]>0) {
        res_de_dw_id_l <- cbind(res_de_dw_id, paste(groupA, "_lowerThan_",groupB, sep="."))
        write.table(as.data.frame(res_de_dw_id_l),
                    file=paste0(output_prefix,".DESeq2.all.DE"),
                    sep="\t",quote=F, row.names=F, col.names=F, append=T)
    }
    
    res_output$level <- ifelse(res_output$padj<=padj,
                               ifelse(res_output$log2FoldChange>=log2FC,
                                      paste(groupA,"UP"),
                                      ifelse(res_output$log2FoldChange<=(-1)*(log2FC),
                                             paste(groupB,"UP"), "NoDiff")) , "NoDiff")
    
    volcanoPlot(res_output, "log2FoldChange", "padj",
                "level", saveplot=paste0(file_base1,".volcano.pdf"), ...)
    
    rankPlot(res_output, label=10, saveplot=paste0(file_base1,".rankplot.pdf"), width=20, ...)
    
    
    res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
    res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))
    
    res_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)
    
    
    res_de_top20_expr <- normalized_counts[res_de_top20,]
    
    sample = as.data.frame(dds@colData)
    
    pheatmap::pheatmap(res_de_top20_expr, cluster_row=T, scale="row",
                       annotation_col=sample,
                       filename=paste0(file_base1,".top20DEgenes.heatmap.pdf"))
    
    res_de_top20_expr2 <- data.frame(Gene=rownames(res_de_top20_expr), res_de_top20_expr)
    res_de_top20_expr2 <- reshape2::melt(res_de_top20_expr, id=c("Gene"))
    
    colnames(res_de_top20_expr2) <- c("Gene", "Sample", "Expression")
    
    res_de_top20_expr2$Group <- sample[match(res_de_top20_expr2$Sample, rownames(sample)),design]
    
    p = ggplot(res_de_top20_expr2, aes(x=Gene, y=Expression)) +
        geom_point(aes(color=Group), alpha=0.5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              axis.title.x = element_blank()) +
        ylab("Normalized xpression value") + scale_y_log10()
    ggsave(p, file=paste0(file_base1,".top20DEgenes.dotplot.pdf"), width=20,
           height=14, units="cm", ...)
    
}





dds <- Pett_readscount2deseq(count_matrix_file = exprSet, sampleFile = group_df, design = "group", covariate=NULL )
is_dds_rlog = T
if(!is_dds_rlog){is_dds_vst=T}else{is_dds_vst=F}
normexpr <- Pett_deseq2normalizedExpr(dds, rlog=is_dds_rlog, vst=is_dds_vst)### 数据集小于30 -> rlog，大数据集 -> VST

if(is_dds_rlog){
    exp_miRNA = normexpr$rlog
}else{
    exp_miRNA = normexpr$vst
}

DEG_miRNA_df = as.data.frame(  results(  dds  )  )

sequencing_dds_normexpr_exp_group_pheno_list = list(dds=dds, DEG_miRNA_df=DEG_miRNA_df, normexpr=normexpr, exp=exp, group=group_df, pheno=pheno )
save(sequencing_dds_normexpr_exp_group_pheno_list, file='./data/miRNA_sequencing_dds_normexpr_exp_group_pheno_list.rdata' )

















