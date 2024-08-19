library(data.table)
library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)

#Read in dfs and TADs
df <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x_0011_df_phyloP_TADs.csv")
unique_TADs <- fread("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

#Remove nested TADs
df.no_overlap<-df[ df$TAD_ID %in% unique_TADs$TAD_ID, ]
df.no_overlap$CHRPOS<-paste0(df.no_overlap$seqnames,":",df.no_overlap$start,"-",df.no_overlap$end)

elements<-c("Candidate Weak Enhancer","Candidate Strong Enhancer","Active Promoter","Inactive Promoter")
element.table<-tibble(Relationship=c("Element in Multiple TADs",
                                     "Element within Gene",#"Element crosses TAD boundary",
                                     "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                     "Element-Gene Within TAD"))

for (el in elements){
  test.element<-df.no_overlap%>%
    filter(is.na(geneStart),chromHMM_cat_longest==el)%>%
    select(c(chromHMM_cat_longest,CHRPOS,start,end,seqnames,TAD_chr,TAD_end,TAD_start,geneStart,geneEnd,TAD_ID,distanceToTSS,geneLength))
  n_occur <- data.frame(table(test.element$CHRPOS))
  test.element<-test.element%>%
    mutate(element_check=case_when(CHRPOS %in% n_occur$Var1[n_occur$Freq > 1]~"Element in Multiple TADs",
                                   start<TAD_start|end>TAD_end~"Element crosses TAD boundary",
                                   TRUE~"No"))%>%
    distinct(CHRPOS,.keep_all=TRUE)%>%
    mutate(gene_check=case_when(geneStart<start&start<geneEnd~"In",
                                geneStart<end&end<geneEnd~"In",
                                TRUE~"Out"),
           Relationship=case_when(element_check=="Element in Multiple TADs"~"Element in Multiple TADs",
                                  gene_check=="In"~"Element within Gene",
                                  #element_check=="Element crosses TAD boundary"~"Element crosses TAD boundary",
                                  geneEnd<TAD_start~"Gene off TAD",
                                  geneStart>TAD_end~"Gene off TAD",
                                  geneStart<TAD_start~"Gene Start off TAD",
                                  geneEnd>TAD_end~"Gene End off TAD",
                                  TRUE~"Element-Gene Within TAD"))
  element.table.hold<-test.element%>%
    group_by(Relationship) %>% 
    summarise(n=n()) %>% 
    mutate(percent=round((n/sum(n)*100),1)) %>%
    mutate(Element=paste0(n," (",percent,"%)"))%>%
    select(Relationship,Element)
  element.table<-merge(element.table,element.table.hold,by=c("Relationship"), all=TRUE)
}

colnames(element.table)<-c("Relationship","Cand. Enhancer (Weak)","Cand. Enhancer (Strong)","Active Promoter","Inactive Promoter")
element.table$Relationship2 <- factor(element.table$Relationship, levels = c("Element within Gene","Element-Gene Within TAD",
                                                                             "Gene off TAD","Gene Start off TAD","Gene End off TAD",
                                                                             "Element in Multiple TADs"))
element.table%>%
  arrange(Relationship2)%>%
  select(-c(Relationship2))%>%
  kable()%>%
  kable_paper() %>% 
  pack_rows(index = c("Proximity Assignement Likely Correct" = 2, "Proximity Assignment Less Likely" = 3, "Proximity Assignment Uncertain" = 1))
  

test.element<-df.no_overlap%>%
  filter(!is.na(geneStart),chromHMM_cat_longest=="Active Promoter")%>%
  select(c(chromHMM_cat_longest,CHRPOS,start,end,seqnames,TAD_chr,TAD_end,TAD_start,geneStart,geneEnd,TAD_ID,distanceToTSS,geneLength))
n_occur <- data.frame(table(test.element$CHRPOS))
test.element<-test.element%>%
  mutate(element_check=case_when(CHRPOS %in% n_occur$Var1[n_occur$Freq > 1]~"Element in Multiple TADs",
                                 start<TAD_start|end>TAD_end~"Element crosses TAD boundary",
                                 TRUE~"No"))%>%
  distinct(CHRPOS,.keep_all=TRUE)%>%
  mutate(gene_check=case_when(geneStart<start&start<geneEnd~"In",
                              geneStart<end&end<geneEnd~"In",
                              TRUE~"Out"),
         Relationship=case_when(element_check=="Element in Multiple TADs"~"Element in Multiple TADs",
                                gene_check=="In"~"Element within Gene",
                                #element_check=="Element crosses TAD boundary"~"Element crosses TAD boundary",
                                geneEnd<TAD_start~"Gene off TAD",
                                geneStart>TAD_end~"Gene off TAD",
                                geneStart<TAD_start~"Gene Start off TAD",
                                geneEnd>TAD_end~"Gene End off TAD",
                                TRUE~"Element-Gene Within TAD"))



--------


import pandas as pd
import numpy as np
from tabulate import tabulate

# Read in dfs and TADs
df = pd.read_csv("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x_0011_df_phyloP_TADs.csv")
unique_TADs = pd.read_csv("C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/TAD_no_nested.csv")

# Remove nested TADs
df_no_overlap = df[df['TAD_ID'].isin(unique_TADs['TAD_ID'])]

#creating Chromosome position again
df_no_overlap['CHRPOS'] = df_no_overlap['seqnames'].astype(str) + ":" + df_no_overlap['start'].astype(str) + "-" + df_no_overlap['end'].astype(str)


#declaring elements
elements = ["Candidate Weak Enhancer", "Candidate Strong Enhancer", "Active Promoter", "Inactive Promoter"]


#declaring relationships 
relationship_categories = ["Element in Multiple TADs", "Element within Gene", 
                           "Gene off TAD", "Gene Start off TAD", "Gene End off TAD", "Element-Gene Within TAD"]

#creating df from relationship categories
element_table = pd.DataFrame({'Relationship': relationship_categories})

#iterating thorugh possible relationships listed above not DHSs!
for el in elements:
    #all elements where geneStart exists and the chromHMM annotation is equal to current element
    test_element = df_no_overlap[(~df_no_overlap['geneStart'].isna()) & (df_no_overlap['chromHMM_cat_longest'] == el)]
    #filtering columns
    test_element = test_element[['chromHMM_cat_longest', 'CHRPOS', 'start', 'end', 'seqnames', 'TAD_chr', 'TAD_end', 'TAD_start', 'geneStart', 'geneEnd', 'TAD_ID', 'distanceToTSS', 'geneLength']]
    

    #code for checking the frequency of chromosome appearance 
    n_occur = test_element['CHRPOS'].value_counts().reset_index()
    n_occur.columns = ['CHRPOS', 'Freq']
    
    #uses the frequency check to see if element in multiple TADS or just overlapping with singular TAD
    test_element['element_check'] = np.where(test_element['CHRPOS'].isin(n_occur['CHRPOS'][n_occur['Freq'] > 1]), "Element in Multiple TADs",
                                             np.where((test_element['start'] < test_element['TAD_start']) | (test_element['end'] > test_element['TAD_end']), "Element crosses TAD boundary", "No"))
    
    #removes duplicates if ther are any (will only occur if  "Element in Multiple TADs"
    test_element = test_element.drop_duplicates(subset='CHRPOS')
    

    #does similar thing for gene where checks if gene within or out of TAD
    test_element['gene_check'] = np.where((test_element['geneStart'] < test_element['start']) & (test_element['start'] < test_element['geneEnd']), "In",
                                          np.where((test_element['geneStart'] < test_element['end']) & (test_element['end'] < test_element['geneEnd']), "In", "Out"))
    
    #conditional here where element relationship decided with gene 
    test_element['Relationship'] = np.select(
        [
            test_element['element_check'] == "Element in Multiple TADs",
            test_element['gene_check'] == "In",
            test_element['geneEnd'] < test_element['TAD_start'],
            test_element['geneStart'] > test_element['TAD_end'],
            test_element['geneStart'] < test_element['TAD_start'],
            test_element['geneEnd'] > test_element['TAD_end']
        ],
        [
            "Element in Multiple TADs",
            "Element within Gene",
            "Gene off TAD",
            "Gene off TAD",
            "Gene Start off TAD",
            "Gene End off TAD"
        ],
        default="Element-Gene Within TAD"
    )
    
    #percentile statistics stuff related to what percentage any particular relationship has 
    element_table_hold = test_element.groupby('Relationship').size().reset_index(name='n')
    element_table_hold['percent'] = round((element_table_hold['n'] / element_table_hold['n'].sum()) * 100, 1)
    element_table_hold['Element'] = element_table_hold['n'].astype(str) + " (" + element_table_hold['percent'].astype(str) + "%)"
    element_table_hold = element_table_hold[['Relationship', 'Element']]
    
    #adds all the relationships to table for this given element type
    element_table = element_table.merge(element_table_hold, on='Relationship', how='outer')

element_table.columns = ["Relationship", "Cand. Enhancer (Weak)", "Cand. Enhancer (Strong)", "Active Promoter", "Inactive Promoter"]

# The following part for kable is not directly translatable to Python as kable is specific to R
# Instead, we will use the tabulate library for a simple table output
print(tabulate(element_table, headers='keys', tablefmt='psql'))