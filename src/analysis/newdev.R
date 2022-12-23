library(tidyr)
library(cowplot)

premutate.PC <- function(dataset) {
	dataset$split="pc"
	dataset$insurance_status='not_reported'
	
	# Translate race to common format
	races_initial <- c('white', 'asian', 'black', 'declined', 'amerian ind', 'other', 'pac islander', 'two races')
	races_final <- c('White', 'Asian', 'Black', 'not_reported', 'Other/multiple', 'Other/multiple', 'Other/multiple', 'Other/multiple')
	
	# Translate ethnicity to common format
	ethnicity_initial <- c("hispanic", "non-hispanic", "declined")
	ethnicity_final <- c("Hispanic or Latino", "Non-Hispanic", "not_reported")
	
	res <- dataset %>%
	  dplyr::mutate(
	    race = races_final[match(dataset$race, races_initial)],
	    ethnicity = ethnicity_final[match(dataset$ethnicity, ethnicity_initial)]
	  )
	
	return(res)
}

premutate.ED <- function(dataset) {
	# Translate race to common format
	races_initial <- c('White or Caucasian', 'Asian', 'Black or African American', 'not_reported', 'American Indian or Alaska Native', 'Hispanic/Latino', 'Native Hawaiian or Other Pacific Islander', 'Other', 'Patient Refused', 'Unknown')
	races_final <- c('White', 'Asian', 'Black', 'not_reported', 'Other/multiple', 'Other/multiple', 'Other/multiple', 'Other/multiple', 'not_reported', 'not_reported')

	# Translate ethnicity to common format
	ethnicity_initial <- c("Hispanic or Latino", "Non-Hispanic", "not_reported", "Patient Refused", "Unknown")
	ethnicity_final <- c("Hispanic or Latino", "Non-Hispanic", "not_reported", "not_reported", "not_reported")
	
	res <- dataset %>%
		dplyr::mutate(
			race = races_final[match(dataset$race, races_initial)],
			ethnicity = ethnicity_final[match(dataset$ethnicity, ethnicity_initial)]
		)
		
	return(res)
}

getDemographics <- function(dataset, premut=NULL) {
	# Specify sort orders
	vkey_order <- c('all', 'UCX_abnormal', 'age', 'gender', 'race', 'ethnicity', 'insured_status')
	vvar_order <- c('all', 'no', 'yes', '[18,25]', '(25,35]', '(35,45]', '(45,55]', '(55,65]', '(65,75]', '(75,999]', 'Male', 'Female', 'Asian', 'Black', 'White', 'Other/multiple', 'Hispanic/Latino', 'not_reported', 'Non-Hispanic', 'Hispanic or Latino', 'Commercial', 'Medicaid', 'Medicare', 'Self pay', 'Other', 'not_reported')

	if( !is.null(premut) ) {
		dataset <- premut(dataset)
	}

	# Get the total number of individuals in the dataset
	res1 <- dataset %>%
		mutate(vvar = 'all', vkey='all') %>%
		dplyr::group_by(split, vkey, vvar) %>%
		dplyr::summarize(count=n()) %>%
		mutate( prty = count ) %>%
		pivot_wider(names_from='split', values_from=c('count', 'prty')) %>%
		dplyr::ungroup()		
	
	# Get the stratification table
	res2 <- dataset %>%
		dplyr::mutate(
			age = cut(dataset$age, breaks=c(18, 25, 35, 45, 55, 65, 75, 999), include.lowest=T)
			#race = races_final[match(dataset$race, races_initial)],
			#ethnicity = ethnicity_final[match(dataset$ethnicity, ethnicity_initial)]
		) %>%
		pivot_longer(c(
				UCX_abnormal,
				age,
				gender,
				race,
				ethnicity,
				insurance_status
		), names_to='vkey', values_to='vvar') %>%
		dplyr::group_by(split, vkey, vvar) %>%
		dplyr::summarize(count=n()) %>%
		dplyr::ungroup() %>%
		dplyr::group_by(split, vkey) %>%
		dplyr::mutate(
			prty = prettyNPct(count, count/sum(count)) 
		) %>%
		pivot_wider(names_from='split', values_from=c('count', 'prty')) %>%
		arrange(
			match(vkey, vkey_order),
			match(vvar, vvar_order)
		) %>%
		
		dplyr::ungroup()
		
	return(rbind(res1, res2))	
}

getStratification <- function(dataset, premut=NULL) {
	# Specify sort orders
	vkey_order <- c(
		'age',
		'gender',
		'dysuria',
		'abd_pain',
		'Urinary_tract_infections',
		'ua_blood',
		'ua_clarity',
		'ua_glucose',
		'ua_ketones',
		'ua_leuk',
		'ua_nitrite',
		'ua_protein',
		'ua_bacteria',
		'ua_epi',
		'ua_wbc'
	)
	vvar_order <- c('1', '0', '[18,25]', '(25,35]', '(35,45]', '(45,55]', '(55,65]', '(65,75]', '(75,999]', 'Male', 'Female', 'negative', 'positive', 'small', 'none', 'few', 'moderate', 'large', '4+', 'clear', 'not_clear', 'Yes', 'No', 'other', 'many', 'marked', 'not_reported')


	if( !is.null(premut) ) {
		dataset <- premut(dataset)
	}

	# Get the stratification table
	res <- dataset %>%
		dplyr::mutate(
			age = cut(dataset$age, breaks=c(18, 25, 35, 45, 55, 65, 75, 999), include.lowest=T)
			#race = races_final[match(dataset$race, races_initial)],
			#ethnicity = ethnicity_final[match(dataset$ethnicity, ethnicity_initial)]
		) %>%
		pivot_longer(c(
				#age,
				#gender,
				#race,
				#ethnicity,
				#insurance_status,				
				age,
				gender,
				dysuria,
				abd_pain,
				Urinary_tract_infections,
				ua_blood,
				ua_clarity,
				ua_glucose,
				ua_ketones,
				ua_leuk,
				ua_nitrite,
				ua_protein,
				ua_bacteria,
				ua_epi,
				ua_wbc
				
		), names_to='vkey', values_to='vvar') %>%
		dplyr::group_by(split, UCX_abnormal, vkey, vvar) %>%
		dplyr::summarize(count=n()) %>%
		dplyr::ungroup() %>%
		dplyr::group_by(split, vkey, vvar) %>%
		dplyr::mutate(
			count = prettyNPct(count, count/sum(count)) 
		) %>%
		ungroup() %>%
		pivot_wider(names_from=c('split', 'UCX_abnormal'), values_from=c('count'), values_fill='-') %>%
		arrange(
			match(vkey, vkey_order),
			match(vvar, vvar_order)
		) %>%
		
		dplyr::ungroup()
		
	return(res)
}

#getEDStats(urine_full)
#getEDStratification(urine_full)

	strat_vkey_order <- c(
		'age',
		'gender',
		'dysuria',
		'abd_pain',
		'Urinary_tract_infections',
		'ua_blood',
		'ua_clarity',
		'ua_glucose',
		'ua_ketones',
		'ua_leuk',
		'ua_nitrite',
		'ua_protein',
		'ua_bacteria',
		'ua_epi',
		'ua_wbc'
	)
	strat_vvar_order <-  c('1', '0', '[18,25]', '(25,35]', '(35,45]', '(45,55]', '(55,65]', '(65,75]', '(75,999]', 'Male', 'Female', 'negative', 'positive', 'small', 'none', 'few', 'moderate', 'large', '4+', 'clear', 'not_clear', 'Yes', 'No', 'other', 'many', 'marked', 'not_reported')


stratTable <- full_join(
		getStratification(urine_full, premut=premutate.ED), 
		getStratification(primarycare_dat, premut=premutate.PC)) %>%
	replace_na(list(
		'training_no'='-', 
		'training_yes'='-',
		'validation_no'='-',
		'validation_yes'='-',
		'pc_no'='-',
		'pc_yes'='-'
	)) %>%
	dplyr::arrange(
		match(vkey, strat_vkey_order),
		match(vvar, strat_vvar_order)
	)
write.table(stratTable, 'stratification-table.txt', quote=F, row.names=F, sep="\t")

plot_grid(
	fig.roc.int,
	calibFig.int, 
	fig.roc.ext,
	calibFig.ext,
	labels=c('A', 'B', 'C', 'D'), ncol=2)
	
ggsave('all.png', height=18, width=18, unit="in", dpi=1200)

demoED <- getDemographics(urine_full, premut=premutate.ED)
demoPC <- getDemographics(primarycare_dat, premut=premutate.PC)

table1 <- full_join(demoED, demoPC)
write.table(table1, 'table1.txt', sep="\t", quote=F, row.names=F)
