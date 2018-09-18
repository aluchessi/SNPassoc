

```R
getwd()
```


'/home/adluchessi/Script'



```R
library(SNPassoc)
```

    Loading required package: haplo.stats
    Loading required package: survival
    Loading required package: mvtnorm
    Loading required package: parallel



```R
db1<- read.table("bd_raul.txt", sep="", header = T)
```


```R
summary(db1)
```


           ID             SEX          AGE        BLOOD     HAIR         Group   
     Min.   :  1.00   Female: 8   Min.   :14.00   No :52   No :50   negative:40  
     1st Qu.: 28.25   male  :98   1st Qu.:23.00   Yes:54   yes:56   Positive:66  
     Median : 54.50               Median :27.00                                  
     Mean   : 54.49               Mean   :31.67                                  
     3rd Qu.: 80.75               3rd Qu.:37.75                                  
     Max.   :107.00               Max.   :68.00                                  
                                                                                 
          Ethinia       Height         Wheight      rs6871510 rs6454674 rs684513 
     non-white:33   Min.   :1.600   Min.   :60.00   CC  :57   GG  : 9   CC  :47  
     white    :25   1st Qu.:1.650   1st Qu.:70.00   CT  :39   GT  :45   CG  :42  
     NA's     :48   Median :1.700   Median :70.00   TT  : 5   TT  :47   GG  : 7  
                    Mean   :1.711   Mean   :73.21   NA's: 5   NA's: 5   NA's:10  
                    3rd Qu.:1.750   3rd Qu.:80.00                                
                    Max.   :1.850   Max.   :95.00                                
                    NA's   :17      NA's   :64                                   



```R
a<-setupSNP(db1,colSNPs=10:12,sep="")
summary(a)
tableHWE(a,Group)

```

              alleles major.allele.freq HWE      missing (%)
    rs6871510 C/T     75.7              0.788033 4.7        
    rs6454674 T/G     68.8              0.818624 4.7        
    rs684513  C/G     70.8              0.804337 9.4        



<table>
<thead><tr><th></th><th scope=col>all groups</th><th scope=col>negative</th><th scope=col>Positive</th></tr></thead>
<tbody>
	<tr><th scope=row>rs6871510</th><td>0.7880331</td><td>0.6612757</td><td>0.4852654</td></tr>
	<tr><th scope=row>rs6454674</th><td>0.8186243</td><td>1.0000000</td><td>0.7761170</td></tr>
	<tr><th scope=row>rs684513</th><td>0.8043374</td><td>1.0000000</td><td>0.7363591</td></tr>
</tbody>
</table>




```R
WGassociation(Group~1,a)
ans<-WGassociation(Group,a)
ans
plot(ans)

```


<table>
<thead><tr><th></th><th scope=col>comments</th><th scope=col>codominant</th><th scope=col>dominant</th><th scope=col>recessive</th><th scope=col>overdominant</th><th scope=col>log-additive</th></tr></thead>
<tbody>
	<tr><th scope=row>rs6871510</th><td>NA       </td><td>0.4792525</td><td>0.6828986</td><td>0.3218007</td><td>0.3853745</td><td>0.9781358</td></tr>
	<tr><th scope=row>rs6454674</th><td>NA       </td><td>0.5052434</td><td>0.2425986</td><td>0.7309163</td><td>0.3272147</td><td>0.2867549</td></tr>
	<tr><th scope=row>rs684513</th><td>NA       </td><td>0.4814842</td><td>0.3744978</td><td>0.3013049</td><td>0.7313941</td><td>0.2517408</td></tr>
</tbody>
</table>




<table>
<thead><tr><th></th><th scope=col>comments</th><th scope=col>codominant</th><th scope=col>dominant</th><th scope=col>recessive</th><th scope=col>overdominant</th><th scope=col>log-additive</th></tr></thead>
<tbody>
	<tr><th scope=row>rs6871510</th><td>NA       </td><td>0.4792525</td><td>0.6828986</td><td>0.3218007</td><td>0.3853745</td><td>0.9781358</td></tr>
	<tr><th scope=row>rs6454674</th><td>NA       </td><td>0.5052434</td><td>0.2425986</td><td>0.7309163</td><td>0.3272147</td><td>0.2867549</td></tr>
	<tr><th scope=row>rs684513</th><td>NA       </td><td>0.4814842</td><td>0.3744978</td><td>0.3013049</td><td>0.7313941</td><td>0.2517408</td></tr>
</tbody>
</table>



    Warning: No SNP is statistically significant after 
             Bonferroni Correction under codominant model 
    Warning: No SNP is statistically significant after 
             Bonferroni Correction under dominant model 
    Warning: No SNP is statistically significant after 
             Bonferroni Correction under recessive model 
    Warning: No SNP is statistically significant after 
             Bonferroni Correction under overdominant model 
    Warning: No SNP is statistically significant after 
             Bonferroni Correction under log-additive model 



![png](output_5_3.png)



```R
association(Group~rs6454674,a)
```


    
    SNP: rs6871510  adjusted by: 
                 negative    % Positive    %   OR lower upper p-value   AIC
    Codominant                                                             
    C/C                23 59.0       34 54.8 1.00              0.4793 139.3
    C/T                13 33.3       26 41.9 1.35  0.58  3.17              
    T/T                 3  7.7        2  3.2 0.45  0.07  2.91              
    Dominant                                                               
    C/C                23 59.0       34 54.8 1.00              0.6829 138.6
    C/T-T/T            16 41.0       28 45.2 1.18  0.53  2.66              
    Recessive                                                              
    C/C-C/T            36 92.3       60 96.8 1.00              0.3218 137.8
    T/T                 3  7.7        2  3.2 0.40  0.06  2.51              
    Overdominant                                                           
    C/C-T/T            26 66.7       36 58.1 1.00              0.3854 138.0
    C/T                13 33.3       26 41.9 1.44  0.63  3.33              
    log-Additive                                                           
    0,1,2              39 38.6       62 61.4 0.99  0.50  1.95  0.9781 138.7



```R
association(Group~rs684513,a)
```


    
    SNP: rs6454674  adjusted by: 
                 negative    % Positive    %   OR lower upper p-value   AIC
    Codominant                                                             
    T/T                21 53.8       26 41.9 1.00              0.5052 139.4
    G/T                15 38.5       30 48.4 1.62  0.69  3.76              
    G/G                 3  7.7        6  9.7 1.62  0.36  7.24              
    Dominant                                                               
    T/T                21 53.8       26 41.9 1.00              0.2426 137.4
    G/T-G/G            18 46.2       36 58.1 1.62  0.72  3.62              
    Recessive                                                              
    T/T-G/T            36 92.3       56 90.3 1.00              0.7309 138.6
    G/G                 3  7.7        6  9.7 1.29  0.30  5.47              
    Overdominant                                                           
    T/T-G/G            24 61.5       32 51.6 1.00              0.3272 137.8
    G/T                15 38.5       30 48.4 1.50  0.66  3.39              
    log-Additive                                                           
    0,1,2              39 38.6       62 61.4 1.41  0.74  2.68  0.2868 137.6



```R
association(Group~rs684513,a)
```


    
    SNP: rs684513  adjusted by: 
                 negative    % Positive    %   OR lower upper p-value   AIC
    Codominant                                                             
    C/C                16 43.2       31 52.5 1.00              0.4815 132.5
    C/G                17 45.9       25 42.4 0.76  0.32  1.80              
    G/G                 4 10.8        3  5.1 0.39  0.08  1.94              
    Dominant                                                               
    C/C                16 43.2       31 52.5 1.00              0.3745 131.2
    C/G-G/G            21 56.8       28 47.5 0.69  0.30  1.57              
    Recessive                                                              
    C/C-C/G            33 89.2       56 94.9 1.00              0.3013 130.9
    G/G                 4 10.8        3  5.1 0.44  0.09  2.10              
    Overdominant                                                           
    C/C-G/G            20 54.1       34 57.6 1.00              0.7314 131.9
    C/G                17 45.9       25 42.4 0.87  0.38  1.98              
    log-Additive                                                           
    0,1,2              37 38.5       59 61.5 0.68  0.35  1.32  0.2517 130.7

