#' Combine and predict CN state
#'
#' Creates a series of human readable interpretations of the calls from the BAF and coverage based statistical models.  The level of detail provided in the interpretation differs depending on the number of SNPs present in each region, with more SNPs allowing for more detailed interpretation.  For example, when only a small number of SNPs are available, only "amplified", "deleted", and "neutral" are selected from.
#' 
#' @param covp The output of coverageFDR.
#' @param bafp The output of allBAF_FDR.
#' @param FDR_Cut The FDR cut-off used in calling CN states.
#' @param SNP_Cutoff How many SNPs should a region have before we can confidently use the BAF data?
#' @param coverageNormCut Regions that have absolute log coverage (tumour/normal) less than this cutoff are considered to be CN neutral.
#' @param coverageExtremeCut Regions that have absolute log coverage (tumour/normal) greater than this cutoff are considered to have extreme CN.  This is only used to determine regions that are homozygously deleted or duplicated if we have significant BAF data in support of this conclusion already.
#' @return A data.frame containing the calls for each type of data and interpretations of these calls.
inferCN_State = function(covp,bafp,FDR_Cut=.01,SNP_Cutoff=10,coverageNormCut = .06, coverageExtremeCut = 0.3){
  #Merge the two tables
  totp = covp
  totp$coverageImbalanceFDR = totp$FDR
  totp$FDR = NULL
  o = match(as.character(interaction(totp$Group,totp$Region,drop=TRUE)),
            as.character(interaction(bafp$Group,bafp$Region,drop=TRUE))
            )
  totp$allelicImbalanceFDR = bafp[o,'FDR']
  #totp$allelicBalanceFDR = bafp[o,'FDR_null']
  totp$N_SNP = bafp[o,'N']
  totp$epsilon = bafp[o,'epsilon']
  totp$fA = bafp[o,'fA']
  totp$Delta = covp$logFC
  totp$aberrationFDR = bafp[o,'combFDR']
  #Set defaults for missing data, where NA is not appropriate
  totp$N_SNP[is.na(totp$N_SNP)]=0
  totp$aberrationFDR[is.na(totp$aberrationFDR)]=totp$coverageImbalanceFDR[is.na(totp$aberrationFDR)]
  #Call Allelic Balance, Imbalance and Coverage Imbalance
  #totp$allelicBalance = totp$allelicBalanceFDR < FDR_Cut
  totp$allelicImbalance = totp$allelicImbalanceFDR < FDR_Cut
  totp$coverageImbalance = totp$coverageImbalanceFDR < FDR_Cut
  totp$aberration = totp$aberrationFDR < FDR_Cut
  #Interpret the results
  #Three cases, no BAF data, limited BAF, lots of BAF data
  #For the moment the no BAF and limited BAF are treated the same, but break them out in case we want to change this in future
  totp$Interpretation = ifelse(totp$N_SNP==0,
                               #When we have no SNPs
                               ifelse(totp$coverageImbalance,
                                      ifelse(totp$Delta < 0,
                                             'Deletion',
                                             'Amplification'
                                             ),
                                      'Neutral'
                                      ),
                               ifelse(totp$N_SNP > SNP_Cutoff,
                                      #When we have a large number of SNPs
                                      ifelse(totp$allelicImbalance,
                                             ifelse(totp$coverageImbalance,
                                                    ifelse(totp$Delta <0,
                                                           'Heterozygous Deletion',
                                                           'Unbalanced Amplification'
                                                           ),
                                                    ifelse(abs(totp$Delta) < coverageNormCut,
                                                           'CNLOH',
                                                           'Neutral'
                                                           )
                                                    ),
                                             ifelse(totp$coverageImbalance & abs(totp$Delta) > coverageExtremeCut,
                                                    ifelse(totp$Delta < 0,
                                                           'Homozygous Deletion',
                                                           'Duplication'
                                                           ),
                                                    'Neutral'
                                                    )
                                             ),
                                      #When we have a small number of SNPs
                                      ifelse(totp$aberration,
                                             ifelse(totp$Delta < 0,
                                                    'Deletion',
                                                    'Amplification'
                                                    ),
                                             'Neutral'
                                             )
                                      )
                               )
  ##The case where we weight both data sources equally
  #totp$Interpretation = ifelse(is.na(totp$allelicImbalance),
  #                         ifelse(totp$coverageImbalance,
  #                           ifelse(totp$logFC<0,
  #                             'Deletion.  No BAF data to tell what kind.',
  #                             'Amplification.  No BAF data to tell what kind.'
  #                           ),
  #                           ifelse(totp$logFC<0,
  #                             "Normal. May be a deletion if you believe the coverage change is significant.  No BAF data available.",
  #                             "Normal. May be an amplification if you believe the coverage change is significant.  No BAF data available."
  #                           )
  #                         ),
  #                         ifelse(totp$allelicImbalance,
  #                           ifelse(totp$coverageImbalance,
  #                             ifelse(totp$logFC<0,
  #                               "Heterozygous Deletion.",
  #                               "Unbalanced Amplification."
  #                             ),
  #                             ifelse(totp$logFC<0,
  #                               "CNLOH. May be a heterozygous deletion if you believe the coverage change is significant.",
  #                               "CNLOH. May be an amplification if you believe the coverage change is significant."
  #                             )
  #                           ),
  #                           ifelse(totp$allelicBalance,
  #                             ifelse(totp$coverageImbalance,
  #                               ifelse(totp$logFC<0,
  #                                 "Homozygous Deletion.",
  #                                 "Duplication."
  #                               ),
  #                               ifelse(totp$logFC<0,
  #                                 "Normal. May be a homozygous deletion if you believe the coverage change is significant (which is unlikely).",
  #                                 "Normal. May be a duplication if you believe the coverage change is significant (which is unlikely)."
  #                               )
  #                             ),
  #                             ifelse(totp$coverageImbalance,
  #                               ifelse(totp$logFC<0,
  #                                 "Deletion.  BAF data inconclusive whether it's homozygous or heterozygous.",
  #                                 "Amplification.  BAF data inconclusive whether it's a duplication or amplification."
  #                               ),
  #                               ifelse(totp$logFC<0,
  #                                 "Normal. May be a deletion of unclear type if you believe the coverage change is significant.",
  #                                 "Normal. May be an amplification of unclear type if you believe the coverage change is significant."
  #                               )
  #                             )
  #                           )
  #                         )
  #                       )
  #totp$InterpretationBAF = ifelse(is.na(totp$allelicImbalance),
  #                           'Unknown.  No BAF data.',
  #                           ifelse(totp$allelicImbalance,
  #                             ifelse(totp$logFC<0,
  #                               'Heterozygous Deletion.',
  #                               'Unbalanced Amplification.'
  #                             ),
  #                             ifelse(totp$allelicBalance,
  #                               ifelse(totp$coverageImbalance,
  #                                 ifelse(totp$logFC<0,
  #                                   'Homozygous Deletion.',
  #                                   'Duplication.'
  #                                 ),
  #                                 'Normal.'
  #                               ),
  #                               'Normal.  But BAF call is not conclusive either way.'
  #                             )
  #                           )
  #                         )
  #totp$InterpretationCov = ifelse(totp$coverageImbalance,
  #                           ifelse(totp$logFC<0,
  #                             'Deletion.',
  #                             'Amplification.'
  #                           ),
  #                           'Normal.'
  #                         )
  #totp$ShortInt = switch(primacy,
  #                       None = gsub(' ','',gsub('^(.+?)\\..*$','\\1',covp$Interpretation)),
  #                       BAF = gsub(' ','',gsub('^(.+?)\\..*$','\\1',covp$InterpretationBAF)),
  #                       Coverage = gsub(' ','',gsub('^(.+?)\\..*$','\\1',covp$InterpretationCov))
  #                       )
  return(totp)
}
