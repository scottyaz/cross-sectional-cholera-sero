## ## This file includes some core and helper functions used in the analyses
## ## for Estimating Cholera Incidece from Cross-Sectional Serology.
## ## Contact azman@jhu.edu for any issues

## ## this should point to the base path
my_root_path <- "."

reload_source <- function(root_path=my_root_path){

    cat("sourcing packages and code \n")
    #library(packrat) ## just in case this isn't in .Rprofile
    #packrat::on()

    ## first load packages
    library(purrr)
    library(plotly)
    library(plyr)
    library(dplyr)
    library(magrittr)
    library(tidyr)
    library(randomForest)
    library(pander)
    library(partykit) ## conditional Random Forests but not really using anymore
    library(ROCR)
    library(broom)
    library(grid)
    library(stringr)
    library(readr)
    library(readxl)
    library(stringr)
    library(grid)
    library(knitr)
    library(kableExtra) # for nicer tables
    library(cvAUC)
    library(flextable)
    library(officer)
    library(scam)
    library(gridExtra)
    library(RColorBrewer)
    library(pairwiseCI)



    ## setting ggplot theme to be global
    theme_set(theme_bw())

    ## then source code
    source(paste0(my_root_path,"/source/base.R"))
    source(paste0(my_root_path,"/source/util.R"))

}


## ROCs for cross validated results (multiple)
make_rocs <- function(preds,truths,k_folds,ribbon=TRUE,title="",annot_auc=TRUE){
    
    votes <- sapply(preds,function(pred)
        apply(pred$individual,1,function(x) mean(x==1)),simplify=FALSE)

    perf_auc <- sapply(1:k_folds,function(x)
        prediction(votes[[x]],truths[[x]]) %>% performance(.,"auc"))
  
    # get the area under the curve for each fold
    auc <- sapply(perf_auc,function(x) x@y.values[[1]])

    
    perf_roc <- sapply(1:k_folds,function(x) {
        prediction(votes[[x]],truths[[x]]) %>% performance(.,"tpr","fpr")},
        simplify = F)

    df <- data.frame(fold = rep(1:k_folds,sapply(perf_roc,function(x) x@x.values %>% unlist %>% length)),
                     fpr=sapply(perf_roc,function(x) x@x.values[[1]],simplify=FALSE) %>% unlist,
                     tpr=sapply(perf_roc,function(x) x@y.values[[1]],simplify=FALSE) %>% unlist)


    cat(sprintf("AUC Range: %.2f-%.2f \n",range(auc)[1],range(auc)[2]))
    cat(sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
    cat(sprintf("AUC Median and Mean: %.2f, %.2f \n",median(auc),mean(auc)))

    df <- df %>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))
    ci_ribbon <- df %>% group_by(fpr_r) %>% dplyr::summarize(med=median(tpr),low_ci=quantile(tpr,.025),upper_ci=quantile(tpr,.975))

    if(ribbon){

        ggplot(data=ci_ribbon) +
        geom_line(aes(x=fpr_r,y=med),col='royalblue',lwd=1.2)+
        geom_ribbon(data=ci_ribbon,aes(x=fpr_r,ymin=low_ci,ymax=upper_ci),alpha=.5,fill='royalblue') +
        xlab('False Positive Rate') + ylab('True Positive Rate') -> gg

    } else {

        ggplot(data=df,aes(x=fpr,y=tpr)) +
        geom_line(aes(group=fold),color=AddAlpha('#2c7fb8',0.9)) +
        xlab('False Positive Rate') + ylab('True Positive Rate') +
        annotate("text", x = 0.5, y = 0.15, size=2,label = title) -> gg

        if(annot_auc) gg <- gg + annotate("text", x = 0.5, y = 0.1,size=2, label = paste0("AUC Range=",round(min(auc),3),"-",round(max(auc),3)))

    }

    return(list(gg,df,auc))

}


# extract sens and spec from predictions
sens_spec <- function(truth,pred){
  conf_matrix <- table(data.frame(truth=truth,predicted=pred))

  if(nrow(conf_matrix)==2){
    sens <- conf_matrix[2,2]/sum(conf_matrix[2,])
    spec <- conf_matrix[1,1]/sum(conf_matrix[1,])
  } else {
    sens <- diag(conf_matrix)/rowSums(conf_matrix)
    spec <- sapply(1:nrow(conf_matrix),function(x) sum(conf_matrix[,-x])/sum(conf_matrix[-x,]))
  }

  return(data.frame(sens=sens,spec=spec))
}



##' helper function to load data
##' based on (wide) data cleaned in data cleaning RMD file
##' @title
##' @param long data in long (TRUE) or wide (FALSE) format?
##' @param type 'cases' for only cases, 'controls' for only controls and 'both' for all
##' @return tbf_df object
##' @author asa
load_data <- function(long=TRUE,type="both",root_path=my_root_path){

    dat <- read_csv("data/public_bangladesh_data.csv")

    ## if we want to make it long
    if(long){
        dat <- dat %>%
            tidyr::gather(assay,value,vibinab:pgctxb_s,plpsm_s,pmctxb_s) %>%
            dplyr::mutate(value=as.numeric(value))
    }

    ## cases only? controls only?
    if(type=="cases"){
        dat <- dat %>% dplyr::filter(!is.na(cx))
    } else if(type == "controls"){
        dat <- dat %>% dplyr::filter(is.na(cx))
    }

    if (!type %in% c("cases","controls","both")){
        warning("sorry don't know that type, returning data for cases AND controls")
    }

    return(dat)
}


##' loads data for random forest models and cut-point analysis
##'
##' @param include_igm
##' @param extra_vars
##' @param long -
##' @param full - if FALSE will drop observations with any missing values for any specified variable
load_data_for_paper <- function(include_igm=FALSE,
                                extra_vars,
                                long=TRUE,
                                full=FALSE){

    dat_full <- load_data(long=FALSE)

    std_vars <- vars <- c("id","lbday2","lbday","vibinab","vibogaw",
                          "plpsa_s","ptctxa_s","pglsp_s","pgctxb_s",
                          "sex","wt","ht","age","o_group","cx",
                          "inf_10","inf_45","inf_100","inf_200","inf_365")


    if(include_igm) {
        vars <- c(std_vars,"plpsm_s","pmctxb_s")
    }

    if(!missing(extra_vars)){
        vars <- c(vars,extra_vars)
    }

    my_vars <- enquo(vars)

    dat <- dat_full %>% select(!!my_vars)

    if(full){
        rc <- dat_full
    } else {
        rc <- dat
        ## only returning complete cases based on std vars
        rc <- rc[complete.cases(rc %>% select(!!enquo(std_vars))),]
    }

    if(long && !full){

        assays <- c('vibinab','vibogaw',vars[which(endsWith(vars,"_s"))])

        assays_q <- quo(assays)
        rc <- rc %>%
        tidyr::gather(assay,value,!!assays_q) %>%
        dplyr::mutate(value=as.numeric(value))

    }

    if (long && full){

        ## grabbing all potential titers
        assays <- c('vibinab','vibogaw',colnames(rc)[which(endsWith(colnames(rc),"_s"))])

        assays_q <- quo(assays)

        rc <- rc %>%
        tidyr::gather(assay,value,!!assays_q) %>%
        dplyr::mutate(value=as.numeric(value))

    }

    return(rc)
}



##' make a reduced dataset with only variables of interest
##' this is a function that really isn't needed anymore
##' but used in some code that I don't feel like breaking
##' @title
##' @param vars_of_interest
##' @param serotypes vector of serotypes to include for cases (1 = ogawa, 2=inaba, 9=controls)
##' @param long
##' @param root_path
##' @return
##' @author
make_dat_reduced <- function(vars_of_interest = c('id','lbday','lbday2',
                                                  'age','sex','o_group',
                                                  'ht','wt',
                                                  'vibinab','vibogaw',
                                                  'plpsa_s','ptctxa_s',
                                                  'pglsp_s','pgctxb_s','cx'),
                             serotypes=c(1,2,9),
                             long=TRUE,
                             root_path=my_root_path){

    my_dat <- load_data_for_paper(long=FALSE) %>%
    dplyr::filter(cx %in% serotypes)

    my_dat <- my_dat[,vars_of_interest]

    if(long){
        my_dat <- my_dat %>% gather(assay,value,vibinab:pgctxb_s)
        my_dat <- my_dat %>% dplyr::mutate(value=as.numeric(value))
    }

    return(my_dat)
}

##' Makes a reduced dataset of cases and returns
##'
##' @title
##' @param splines TRUE if we want to sue interpolation splines for trejctories
##' @param ... options for make_dat_reduced_cases
##' @return a list with interpoliation functions for each assay (and a vector called ids with the study ids for each)
##' @author asa
make_all_trajs <- function(serotypes=c(1,2),...){

    dat_l <- make_dat_reduced(serotypes=serotypes,...)

    ## removing those cases with only one observation
    ## which all happen to be on day 2-5 anyway
    rems <- dat_l %>%
    group_by(id,assay) %>%
    summarize(remove=n()<2) %>%
    filter(remove==TRUE) %>%
    select(id) %>% distinct %>% unlist

    dat_l <- dat_l %>% filter(!id %in% rems)

    ## figure out the maximum day observation for each person
    ## to use as a lookup
    look_up_tab <- dat_l %>%
    group_by(id) %>%
    summarize(min_day = min(lbday2),
              max_day = max(lbday2)
              )

    dat_l <- left_join(dat_l,look_up_tab)


    ## now we generate linear interpolations
    ## and splines from data
    ## and create a tibble with these functions embedded (along with the demographic data)
    rc <- dat_l %>%
    group_by(id,assay,min_day,max_day,age,sex,o_group,ht,wt,cx) %>%
    do(linear_func = approxfun(x=.$lbday2,y=.$value),
       spline_func = splinefun(x=.$lbday2, y=.$value,method='monoH.FC')
       )

    return(rc)

}


##' @param t vector with times since infection
##' @param traj_tbl tbl as generated from xxx
##' @param round_vibs will round decimal titers to nearest potential whole titer value
##' @return
##' @author Andrew A
generate_seroprof_dataset <- function(t,traj_tbl,round_vibs=TRUE){

    ## seperating times into Inf, which we will get from
    ## controls and finite t's, which we wil get from case
    ## trajectories
    infs_exist <- any(is.infinite(t))

    if(infs_exist){

        n_t_infs <- sum(is.infinite(t))

        ## generate uninfected profiles for these folks
        dat_uninfs <- generate_uninfected_profile(n_t_infs)
        dat_uninfs$t_since_inf <- Inf

        t <- t[which(is.finite(t))]
    }

    ## if only Infs then we are done
    if(length(t)==0){
        return(dat_uninfs)
    }

    ## check that infection times are actually within the bounds of an
    ## observation

    ## this should be called not_in_obs_range but not changing now
    in_obs_range <- (t<min(traj_tbl$min_day)|t>max(traj_tbl$max_day))
    if(any(in_obs_range)){
        t <- t[!in_obs_range]
        warning("at least one t is not within the range of observed times. dropping these times.")
    }

    ## sample based on available traj data for each person
    ids <- sapply(t,function(inf_time) {
        traj_tbl %>%
        filter(min_day<=inf_time & max_day>=inf_time & assay=="vibogaw") %>%
        sample_n(size=1) %>%
        select(id) %>% unname %>% unlist
    })

    rc <- lapply(1:length(ids),function(x){
        ## supressing warnings about rowwise operations as
        ## they aren't meaningful here
        tmp <- suppressWarnings(
            traj_tbl %>%
            filter(id==ids[x]) %>%
            group_by(assay) %>%
            summarize(value=linear_func[[1]](t[x]))
        )

        dems <- traj_tbl %>%
        filter(id==ids[x]) %>%
        select(id,age:wt) %>%
        slice(1)

        dems$t_since_inf <- t[x]

        bind_cols(dems,
                  tmp %>%
                  select(value) %>% t %>% data.frame %>%
                  set_names(tmp$assay))

    }) %>% bind_rows

    ## add in the uninfected data
    if(infs_exist){
        rc <- bind_rows(rc,dat_uninfs)
    }


    ## for now rounding up to the nearest proper titer
    if(round_vibs){
        potential_vibs <- c(5,10,10*2^(1:20))
        rc <- rc %>%
        rowwise %>%
        mutate(vibogaw = potential_vibs[which(potential_vibs >= vibogaw)[1]],
               vibinab = potential_vibs[which(potential_vibs >= vibinab)[1]])
    }

    return(rc)

}

## generates uninfected demographic and serologic
## profile for n individuals
generate_uninfected_profile <- function(n){

    contacts <- make_dat_reduced(serotypes = c(9),long=F)
    early_and_late_measures <- make_dat_reduced(serotypes = c(1,2),long=F) %>%
    filter(lbday2<5|lbday2>365)
    rc <- bind_rows(contacts,early_and_late_measures)

    rc <- rc[sample(nrow(rc),n,replace=T),]

    rc <- rc %>% rename(t_since_inf=lbday2) %>%
    select(id,age,sex,o_group,ht,wt,t_since_inf,
           pgctxb_s,pglsp_s,plpsa_s,ptctxa_s,
           vibinab,vibogaw)

    return(rc)
}

##' gets single biomarker sens and spec
##' @param dat
##' @param bm
##' @param pv - TRUE if using pax vac data
##' @return
##' @author Andrew A
get_single_bm_sens_spec <- function(dat,bm,pv=FALSE){

    if(pv){
        times <- c('inf_45','inf_100','inf_200')
    } else {
        times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')
    }

    ## when we are using IgM BMs we have some missing data
    bmq <- rlang::sym(bm) #enquo(bm)
    dat <- dat %>% filter(!is.na(!!bmq))

    ## first get cut-offs
    cos <- sapply(times,function(x) test_thresholds(dat,
                                                    biomarker=bm,
                                                    my_threshold = NA,
                                                    col_of_interest=x,
                                                    nsims=1000,
                                                    print_me=FALSE,
                                                    training_prop=0.7) %>% mutate(time=x),
                  simplify=FALSE) %>%
    do.call('rbind',.)

    modal_titers <- cos %>%
    mutate(time=ordered(time,levels= times)) %>%
    group_by(time) %>%
    dplyr::summarize(modal_titer=names(table(threshold))[which.max(table(threshold))] %>%
                     as.numeric)

    ## now use these modal values for performance
    co_perf <- sapply(seq_along(times),
                      function(x) test_thresholds(dat,
                                                  biomarker=bm,
                                                  my_threshold = modal_titers[x,2] %>% unlist,
                                                  col_of_interest=times[x],
                                                  nsims=1000,
                                                  print_me=FALSE,
                                                  training_prop=0.7) %>%
                                  mutate(time=times[x],antibody=bm),simplify=FALSE) %>%
    bind_rows

    ## and make summary table
    sum_tab <- co_perf %>%
    mutate(time=ordered(time,levels=times)) %>%
    filter(!is.nan(value)) %>%
    group_by(time,variable) %>%
    dplyr::summarize(median=median(value),
                     mean=mean(value),
                     sd=sd(value),
                     p2.5=quantile(value,c(.025)),
                     p97.5=quantile(value,c(0.975))) %>%
    left_join(.,modal_titers)


    return(list(summary=sum_tab,full_dist=co_perf))
}

##' adds binary columns for infected in the last n days
##' based on lbday2
##' @title
##' @param dat
##' @return
##' @author asa
add_binary_inf_columns <- function(dat){
    dat %>% dplyr::mutate(
                       inf_10=ifelse(lbday2>4 & lbday<=10 & is_case==1,1,0),
                       inf_45=ifelse(lbday2>4 & lbday2<=45 & is_case==1,1,0),
                       inf_100=ifelse(lbday2>4 & lbday2<=100 & is_case==1,1,0),
                       inf_200=ifelse(lbday2>4 & lbday2<=200 & is_case==1,1,0),
                       inf_365=ifelse(lbday2>4 & lbday2<=365 & is_case==1,1,0))
}

## helper lookup
assay_names_labels <- assay_names <- c(`plpsa_s` = "Anti-LPS IgA",
                        `ptctxa_s` = "Anti-CTB IgA",
                        `vibinab` = "Vibriocidal (Inaba)",
                        `vibogaw` = "Vibriocidal (Ogawa)",
                        `pglsp_s` = "Anti-LPS IgG",
                        `pgctxb_s` = "Anti-CTB IgG",
                        `plpsm_s` = "Anti-LPS IgM",
                        `pmctxb_s` = "Anti-CTB IgM",
                        `vibo_to_lps`= "Vibriocidal to anti-CTB IgG",
                        `vibo_to_ctx` = "Vibriocidal to anti-LPS IgG",
                        `2` = "Day 2",
                        `7` = "Day 7",
                        `30` = "Day 30",
                        `90` = "Day 90",
                        `180` = "Day 180",
                        `360` = "Day 365",
                        `365` = "Day 365",
                        `540` = "Day 540",
                        `720` = "Day 720",
                        `900` = "Day 900",
                        `10` = "Day 10",
                        `28` = "Day 28",
                        `170` = "Day 170",
                        `0` = "Day 0"
                        )

##' Very specific helper function for making ROC curve and variable
##' importance plot from simulations
##' @title
##' @param preds predictions from cross validation
##' @param truths truths (binary)
##' @param imps importance objects
##' @param my_title plot title
##' @return multiplot object
##' @author
make_roc_varimp_plot <- function(preds,truths,imps,my_title,panel_label="",sub=FALSE,ribbon=TRUE,...){

    my_roc <- make_rocs(preds,truths,length(preds),title=my_title,ribbon=ribbon,...)

    importance_df <- do.call('rbind',imps) %>% as_tibble()
    varnames <- rownames(do.call('rbind',imps)) %>% pretty_antibody()
    importance_df$variable <- varnames
    
    importance_df <- importance_df %>% group_by(variable) %>% summarise(MeanDecreaseAccuracy=median(MeanDecreaseAccuracy))

    levels_vars <- importance_df %>%
                   group_by(variable) %>%
                   dplyr::summarise(med=median(MeanDecreaseAccuracy)) %>%
                   select(med) %>% unlist() %>% order()

    impplot <- importance_df %>%
    dplyr::mutate(variable1=factor(variable,
                                   levels=sort(unique(importance_df$variable))[levels_vars])) %>%
    geom_bar(aes(x=variable1,y=MeanDecreaseAccuracy,fill=variable),stat='identity',alpha=0.5) +
    coord_flip() + scale_fill_discrete(guide=F) + xlab('') + ylab('relative importance') +
    theme_classic() + theme(axis.title=element_text(size=8),
                            axis.text.x = element_text(size=6))

    if(sub){
        gg <- my_roc[[1]] + annotate("text",x=0,y=0,label=my_title,hjust=0) +
        annotate("text",x=0,y=1,label=panel_label,hjust=0) +
        annotation_custom(ggplotGrob(impplot), xmin=0.165, xmax=1, ymin=.1, ymax=.7)

        return(gg)

    } else{

        return(multiplot(impplot,my_roc[[1]],cols=2))

    }
}

## some helper functions
sens_func <- function(x,truths){
    mean((truths==x)[which(truths==1)])
}

spec_func <- function(x,truths){
    mean((truths==x)[which(truths==0)])
}

youden_func <- function(x,truths,w_spec=.5){
    (1-w_spec)*sens_func(x,truths) + w_spec*spec_func(x,truths) - 1
}


##' Uses a cross validation approach to get the opimimal cutpoints
##' then tests performance on hold-out. Only uses one observation per
##' person
##' @param dat
##' @param my_threshold threshold for postivie/negative (if known), otehrwise will compute based on weighted youden index
##' @param w_spec weight for youden index
##' @param biomarker
##' @param col_of_interest
##' @param nsims
##' @param print_me
##' @param training_prop
##' @return
##' @author
test_thresholds <- function(dat,
                            my_threshold=NA,
                            w_spec=.5,
                            biomarker='vib',
                            col_of_interest='inf_100',
                            nsims=1000,
                            print_me=FALSE,
                            training_prop=0.7){


    if(biomarker=='vib'){
        ## get unique values of the titer
        unique_titers <- c(dat$vibinab,dat$vibogaw) %>% unique %>% sort

        ## add new column to dat for max_vib
        dat <- dat %>% dplyr::mutate(vib=pmax(vibogaw,vibinab))

    } else {
        if(!biomarker %in% colnames(dat)) stop("biomarker is not a valid name of a column in dat. Only non-column name allowed is 'vib,' which is for max vibriocidal")

        unique_titers <- round(dat[,biomarker]) %>% unlist %>% unique %>% sort
    }

    unique_ids <- unique(dat$id) ## unique ids in the dataset
    obs_indices <- sapply(1:length(unique_ids),
                          function(x) which(dat$id == unique_ids[x])) ## list of observations per person

    ## number of people in training set
    n_train <- (length(unique_ids) * training_prop) %>% ceiling

    truths <- preds <- vector(mode='list', length=nsims)
    testing_perf <- matrix(ncol=5,nrow=nsims) %>% data.frame()
    colnames(testing_perf) <- c('threshold','sens_test','spec_test','sens_train','spec_train')

    for (i in 1:nsims){

        train_ids <- sample(length(unique_ids),n_train,replace=F) ## get ids of those who we want to sample for training
        test_ids <- setdiff(1:length(unique_ids),train_ids) ## then those for testing

        ## sample one time point for each person
        train_obs_inds <- sapply(train_ids,function(x) sample(obs_indices[[x]],1))
        test_obs_inds <- sapply(test_ids,function(x) sample(obs_indices[[x]],1))

        training_titers <- dat[train_obs_inds,biomarker]
        training_truths <- dat[train_obs_inds,col_of_interest]

        testing_titers <- dat[test_obs_inds,biomarker]
        testing_truths <- dat[test_obs_inds,col_of_interest]

        ## now get optimal cut-point (if needed) from the
        ## training set
        if(is.na(my_threshold)){
            ## returns a matrix with columns for each possible cutpoint and rows for observations
            preds_training <- sapply(unique_titers,function(x) training_titers>=x)

            youden_max_id <-  which.max(apply(preds_training,2,function(x)
                youden_func(x,w_spec=w_spec,truths=training_truths)))

            sens_spec_thresh <- c(unique_titers[youden_max_id],
                                  apply(preds_training,2,function(x) sens_func(x,truths=training_truths))[youden_max_id],
                                  apply(preds_training,2,function(x) spec_func(x,truths=training_truths))[youden_max_id])

            my_thresh <- sens_spec_thresh[1]
            sens_train <- sens_spec_thresh[2]
            spec_train <- sens_spec_thresh[3]

        } else {

            my_thresh <- my_threshold
            preds_training <- training_titers>=my_thresh

            sens_train <- sens_func(preds_training,truths=training_truths)
            spec_train <- spec_func(preds_training,truths=training_truths)

        }

        ## now using the testing data, what is the sens and spec from these cut-offs?
        testing_preds <- data.frame(pred=testing_titers >= my_thresh,
                                    truth=testing_truths)

        testing_perf[i,] <- c(threshold=my_thresh,
                              sens_test=sens_func(testing_preds[,1],testing_preds[,2]),
                              spec_test=spec_func(testing_preds[,1],testing_preds[,2]),
                              sens_train=sens_train,
                              spec_train=spec_train)
    }

    ## don't want to keep this loaded
    library(reshape2)
    thresh_df <- melt(testing_perf,id.vars='threshold')
    detach('package:reshape2')

    if (print_me){
        thresh_df %>% ggplot() +
        geom_histogram(aes(value,fill=factor(threshold)),position='identity') +
        facet_wrap(~variable) +
        labs(title=col_of_interest)  -> gg
        print(gg)
    }

    return(thresh_df)
}


##' one epidemic reconstruction from serosurvey
##' this one uses a single model which you need to feed to it
##' @param n sample size
##' @param times_since_infection_vec times since infection (if n is not provided, we will assume we are interested in all of tese people)
##' @param traj_tbl
##' @param my_forest fitted model
##' @param vib_thresh
##' @return
##' @author asa
one_sim_singlemodel_v2 <- function(n,
                                   pop_profile,
                                   my_forest,
                                   vib_thresh=320,
                                   round_vibs=TRUE){


    sampled_people <- sample(1:nrow(pop_profile),size=n)
    epi_trajs <- pop_profile[sampled_people,]

    ## predict at different time frame
    my_preds <- predict(my_forest,newdata=epi_trajs)

    vib_preds <- epi_trajs %>%
    mutate(pos=ifelse(vibinab>=vib_thresh|vibogaw>=vib_thresh,1,0))

    if(any(is.na(vib_preds))) recover()

    ## final size preds
    rc <- data.frame(
        n=n,
        positive_rf=table(my_preds)[2],
        cum_inc_rf=table(my_preds)[2]/sum(table(my_preds)),
        positive_vib=sum(vib_preds$pos),
        cum_inc_vib=sum(vib_preds$pos)/nrow(vib_preds),
        true_positive_N = sum(pop_profile$t_since_inf <= 365),
        true_positive_n = sum(epi_trajs$t_since_inf <= 365),
        stringsAsFactors = FALSE
    )

    return(rc)
}


##' helper function to make pretty names for covariates
##' either give a dataframe so it works on column names, otherwise
##' give a vector and it replaces names lnly
##' @param x
##' @param dplyr::rename_cols
##' @return
##' @author asa
pretty_antibody <- function(x){

    if(!is.null(ncol(x)) && ncol(x)>1){
        x <- x %>%
        dplyr::rename('plpsa_s',"anti-LPS IgA") %>%
        dplyr::rename('ptctxa_s',"anti-CTB IgA") %>%
        dplyr::rename('vibinab' ,"vibriocidal (Inaba)") %>%
        dplyr::rename('vibogaw' ,"vibriocidal (Ogawa)") %>%
        dplyr::rename('pglsp_s' ,"anti-LPS IgG") %>%
        dplyr::rename('pgctxb_s', "anti-CTB IgG") %>%
        dplyr::rename('plpsm_s',"anti-LPS IgM") %>%
        dplyr::rename('pmctxb_s',"anti-CTB IgM") %>%
        dplyr::rename('lbday2',"Days from Symptom Onset") %>%
        dplyr::rename('wt' ,"weight") %>%
        dplyr::rename('ht' ,"height") %>%
        dplyr::rename('o_group' ,"O blood group")

    } else {

    x <- x %>%
    str_replace_all('plpsa_s',"anti-LPS IgA") %>%
    str_replace_all('ptctxa_s',"anti-CTB IgA") %>%
    str_replace_all('vibinab' ,"vibriocidal (Inaba)") %>%
    str_replace_all('vibogaw' ,"vibriocidal (Ogawa)") %>%
    str_replace_all('pglsp_s' ,"anti-LPS IgG") %>%
    str_replace_all('pgctxb_s', "anti-CTB IgG") %>%
    str_replace_all('plpsm_s' ,"anti-LPS IgM") %>%
    str_replace_all('pmctxb_s' ,"anti-CTB IgM") %>%
    str_replace_all('lbday2' ,"Days from Symptom Onset") %>%
    str_replace_all('^wt$' ,"weight") %>%
    str_replace_all('^ht$' ,"height") %>%
    str_replace_all('o_group' ,"O blood group")

    }
    return(x)
}

## creates folds for cross validation
## for repeated measures
## keeps individuals together within folds
## but no balance by outcome
## returns a list with the observation indexs
## for each fold
cv_folds_rm <- function(dat, k){

    ## get number of observations per person
    unique_ids <- unique(dat$id) ## unique ids in the dataset
    obs_indices <- lapply(1:length(unique_ids),function(x) which(dat$id == unique_ids[x])) ## list of observations per person
    ## get warnings for not having equal values in splits
    person_id_folds <- suppressWarnings(split(unique_ids,rep(1:k,length(unique_ids))))
    obs_id_folds <- lapply(person_id_folds,
                           function(my_fold){which(dat$id %in% my_fold)})

    return(obs_id_folds)
}



## cross validate single titer measuresment
cv_titer_repeated_measure_auc <- function(dat,k,biomarker,col_of_interest){
    library(cvAUC)

    titers <- dat[,biomarker] %>% unlist
    truths <- dat[,col_of_interest] %>% unlist
    ids <- dat[,'id'] %>% unlist

    folds <- cv_folds_rm(dat,k)
    preds <- vector("list",length=k)

    for (k in 1:length(folds)){
        perf <- ROCR::performance(ROCR::prediction(titers[-folds[[k]]],truths[-folds[[k]]]), "sens", "spec")
        df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
        cutoff <- df[which.max(df$sens + df$spec), "cut"]
        preds[[k]] <- ifelse(titers[folds[[k]]]>=cutoff,1,0) %>% unlist %>% unname
    }

    out <- ci.pooled.cvAUC(predictions=preds %>% unlist,
                           labels=truths[unlist(folds)],
                           ids=ids[unlist(folds)],
                           folds=rep(1:k,map(folds,length)),
                           confidence=0.95)

    return(out)
}

## cross validated AUC for RF models
cv_rf_model <- function(dat,my_formula,k=20,ntrees=2000,...){

    ## creating folds
    folds <- cv_folds_rm(dat,k)

    my_preds <- my_truths <- numeric()
    window <- str_extract(my_formula[2] %>% as.character,"inf_[0-9_]+")
    tmp_pred <- var_importance <- vector("list",length=k)

    for (f in 1:length(folds)){
        cat(sprintf("Cross-validating AUC, fold %s \n",f))
        my_forest <- randomForest(my_formula,
                                  ntree=ntrees,
                                  data=dat[-c(folds[[f]]),],
                                  importance=TRUE,...)

        my_truths <- c(my_truths,dat[folds[[f]],window] %>% unlist)
        tmp_pred[[f]] <- predict(my_forest,newdata=dat[c(folds[[f]]),],predict.all=TRUE)
        var_importance[[f]] <- importance(my_forest)

        ## get probability from trees
        my_preds <- c(my_preds,
                      apply(tmp_pred[[f]][[2]],1,function(x) mean(x==1)) %>%
                      as.numeric
                      )
    }

    out <- vector("list",length=4)

    out[[1]] <- ci.pooled.cvAUC(predictions=my_preds,
                                labels=my_truths,
                                ids=dat[unlist(folds),'id'] %>% unlist,
                                folds=rep(1:k,map(folds,length)),
                                confidence=0.95) %>%
    unlist %>% t %>%
    data.frame

    out[[1]] <- out[[1]] %>% mutate(time_window=window)

    out[[2]] <- bind_cols(truth=my_truths,pred=my_preds,fold=rep(1:k,map(folds,length)))

    out[[3]] <- tmp_pred

    out[[4]] <- var_importance

    names(out) <- c("perf_summary","outs","preds_full","var_imp")
    return(out)
}

##' correct proportion by known sens and spec
##' note that we return 0 when corrected prop is negative
##' @param num_pos
##' @param num_neg
##' @param sens
##' @param spec
##' @return
##' @author Andrew A
correct_sero_misclass <- function(num_pos,num_neg,sens=.793,spec=.929){
    pmax((num_pos + (num_pos+num_neg)*(spec -1)) / (sens + spec - 1),0)
}

##' helper for scams and plots
##' not very generalizable but used to make figures for paper
##' @title
##' @param assay
##' @param trans_func
##' @param my_igm - including igm?
##' @param ymax
##' @param ymin
##' @param ybreaks
##' @param ylabs
##' @return
##' @author Andrew A
scam_and_plot <- function(assay,
                          trans_func = function(x){log2(x)},
                          my_igm = FALSE,
                          ymax=15,
                          ymin = 0,
                          ybreaks=2:14,
                          ylabs=2^(2:14)){

    my_assay <- enquo(assay)

    dat <- load_data_for_paper(long=FALSE,include_igm = my_igm) %>% ungroup %>%
    filter(!grepl("\\.[0-9]+$",id)) %>%
    mutate(u5=ifelse(age<5,"under 5 years","5 years and older") %>% factor)%>%
    select(id,lbday2,!!my_assay,age,u5) %>%
    rename(titer = !!my_assay) %>%
    mutate(titer=trans_func(titer))

    ## getting data in a similar way to display baseline titers
    dat_baseline <- load_data_for_paper(long=FALSE,include_igm = my_igm) %>% ungroup %>%
    filter(!grepl("\\.[0-9]+$",id)) %>%
    mutate(u5=ifelse(age<5,"under 5 years","5 years and older") %>% factor)%>%
    select(id,lbday,!!my_assay,age,u5) %>%
    rename(titer = !!my_assay) %>%
    mutate(titer=trans_func(titer)) %>%
    filter(lbday==2)

    ## for each person figure out when in the first 30 days their peak
    ## observation is made and then truncate before

    declines <- left_join(dat,
                          dat %>% group_by(id) %>% filter(lbday2<40,!is.na(titer)) %>%
                          summarize(peak_day=max(titer,na.rm=T))
                          ) %>%
    filter(lbday2>=peak_day)


    ## using SCAM models
    scam_dat <- declines
    scam_dat$myid <- scam_dat$id %>% as.factor
    scam_dat <- scam_dat %>% mutate(dum = 1)
    #dum <- rep(1,nrow(scam_dat)) ## dummy variable to help with plotting fixed effects only

    print("running scam model with no covaraiates")

    scam1 <- scam::scam(titer ~ s(lbday2,bs="mpd") +
                  s(myid,bs="re",by=dum),
                  data = scam_dat)

    ## just choosing random person's id to use that RE for plotting
    testdata_1 <- data.frame(lbday2 = seq(2, 900, length = 200),
                             myid="623",
                             dum=0)

    fits_1 <- predict(scam1, newdata=testdata_1, type='response', se=T)
    predicts_1 <- data.frame(testdata_1, fits_1) %>%
    mutate(lower = fit - 1.96*se.fit,
           upper = fit + 1.96*se.fit,
           ymax =quantile(dat_baseline$titer,.75,na.rm=T),
           ymin= quantile(dat_baseline$titer,.25,na.rm=T))


    gg1 <- ggplot(aes(x=lbday2,y=fit), data=predicts_1) +
    geom_ribbon(aes(ymin = lower, ymax=upper), fill='gray90') +
    geom_line() +
    scale_y_continuous(breaks=ybreaks,labels=ylabs) +
    coord_cartesian(ylim=c(ymin, ymax)) +
    theme_light() +
    guides(color=FALSE) +
    xlab("days from symptom onset") +
    ylab("titer") + theme(legend.title=element_blank())

    ## add alpha rectangle with IQR of enrollment titers for cases
    gg1 <- gg1 +
      geom_ribbon(aes(ymax=ymax,ymin=ymin),fill="grey70",alpha=.2)

    print("running scam with covariates (by age)")
    scam2 <- scam(titer ~ s(lbday2,by=u5,bs="mpd") + 
                  s(myid,bs="re",by=dum),
                  family=gaussian(link="identity"),
                  data = scam_dat)

    testdata_2 <- data.frame(lbday2 = rep(seq(2, 900, length = 100),2),
                            u5=rep(c("under 5 years","5 years and older"),each=100),
                            myid = "623",
                            dum=0)

    fits_2 <- predict(scam2, newdata=testdata_2, type='response', se=T)
    predicts_2 = data.frame(testdata_2, fits_2) %>%
    mutate(lower = pmax(0.1,fit - 1.96*se.fit),
           upper = fit + 1.96*se.fit)

    gg2 <- ggplot(data=predicts_2) +
    geom_ribbon(aes(x=lbday2,ymin = lower, ymax=upper,fill=u5),alpha=.5) +
    geom_line(aes(x=lbday2,y=fit,color=u5)) +
    scale_y_continuous(breaks=ybreaks,labels=ylabs) +
    theme_light() +
    coord_cartesian(ylim=c(ymin, ymax)) +
    theme(legend.position = c(.7,.8)) +
    guides(color=FALSE) +
    xlab("days from symptom onset") +
    ylab("titer") + theme(legend.title=element_blank())

    ## need to implement geom_ribbon for each group here
    ribbon_dat <- rbind(
        data.frame(upper=(dat_baseline %>% filter(u5=="under 5 years") %>% summarize(lower=quantile(titer,.75,na.rm=T)) %>% unlist %>% unname),
                   lower=(dat_baseline %>% filter(u5=="under 5 years") %>% summarize(lower=quantile(titer,.25,na.rm=T)) %>% unlist %>% unname),
                   u5 =  "under 5 years",
                   lbday2=2),
        data.frame(upper=(dat_baseline %>% filter(u5=="5 years and older") %>% summarize(lower=quantile(titer,.75,na.rm=T)) %>% unlist %>% unname),
                   lower=(dat_baseline %>% filter(u5=="5 years and older") %>% summarize(lower=quantile(titer,.25,na.rm=T)) %>% unlist %>% unname),
                   u5 = "5 years and older",
                   lbday2=2),
        data.frame(upper=(dat_baseline %>% filter(u5=="under 5 years") %>% summarize(lower=quantile(titer,.75,na.rm=T)) %>% unlist %>% unname),
                   lower=(dat_baseline %>% filter(u5=="under 5 years") %>% summarize(lower=quantile(titer,.25,na.rm=T)) %>% unlist %>% unname),
                   u5 =  "under 5 years",
                   lbday2=900),
        data.frame(upper=(dat_baseline %>% filter(u5=="5 years and older") %>% summarize(lower=quantile(titer,.75,na.rm=T)) %>% unlist %>% unname),
                   lower=(dat_baseline %>% filter(u5=="5 years and older") %>% summarize(lower=quantile(titer,.25,na.rm=T)) %>% unlist %>% unname),
                   u5 = "5 years and older",
                   lbday2=900))

    gg2 <- gg2 + geom_ribbon(data=ribbon_dat,aes(x=lbday2,ymax=upper,ymin=lower,fill=u5),alpha=.15) +
    theme(legend.background=element_blank())

    return(list(gg1,gg2,AIC(scam1,scam2)))
}


##' makes figure 1 with trajectories
##' includes IgM subset
##' @title
##' @param interactive ggplotly version if you want
##' @return
##' @author Andrew A
make_figure_1 <- function(interactive=FALSE){


    dat_l_f <- load_data_for_paper(long=TRUE,full=TRUE,include_igm = TRUE) %>%
    filter(assay %in% c('vibinab','vibogaw',
                        'pglsp_s','pgctxb_s',
                        'plpsa_s','ptctxa_s',
                        'plpsm_s','pmctxb_s')) %>%
    filter(cx!=9) %>%
    mutate(assay = ordered(assay,levels = c('vibinab','vibogaw',
                                    'pglsp_s','pgctxb_s',
                                    'plpsa_s','ptctxa_s',
                                    'plpsm_s','pmctxb_s')))


    ## panel labels
    pls <- data.frame(assay = factor(c('vibinab','vibogaw',
                                       'pglsp_s','pgctxb_s',
                                       'plpsa_s','ptctxa_s',
                                       'plpsm_s','pmctxb_s'),
                                     levels = c('vibinab','vibogaw',
                                                'pglsp_s','pgctxb_s',
                                                'plpsa_s','ptctxa_s',
                                                'plpsm_s','pmctxb_s')),
                      y = rep(Inf,8),
                      x = rep(2,8),
                      hjustvar = rep(0,8) ,
                      vjustvar = rep(1.1,8),
                      mylab = c("A","B","C","D","E","F","G","H"))


    myColors <- brewer.pal(8,"Dark2")
    names(myColors) <- levels(dat_l_f$assay)
    colScale <- scale_colour_manual(name = "assay",values = myColors)

    dat_l_f %>%
    ggplot(aes(x=lbday2,y=value,group=id,color=assay,label=id)) +
    geom_line(alpha=.1) +
    geom_smooth(method='loess',aes(group=assay),se=FALSE) +
    scale_y_continuous(trans='log2') +
    scale_x_log10(breaks=c(2,7,10,30,90,180,365,540,900)) +
    guides(colour=FALSE) + xlab('time since symptom onset (days)') + ylab('titer') +
    geom_text(data=pls,aes(x=x,y=y,label=mylab,hjust=hjustvar,vjust=vjustvar),inherit.aes = FALSE) +
    facet_wrap(~assay,ncol=2,scale="free",labeller = as_labeller(assay_names_labels)) + theme_minimal() +
    colScale-> gg

    if(interactive) {
        return(ggplotly(gg))
    } else {
        return(gg)
    }
}


##' makes figure 2
##' density plots for vibriocidals
##' @return gg object
##' @author Andrew A
make_figure_2 <- function(){

    gg <- make_titer_dist_density_plots(assays=c('vibinab','vibogaw'),
                                        exc_lbdays = c(0,270,720))

    gg <- gg + xlab("vibriocidal titer") +
    scale_x_continuous(trans="log2",breaks= c(5,10,10*2^(1:20))) +
    theme(axis.ticks = element_blank(), legend.background = element_blank(),
          legend.key = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), strip.background = element_blank(),
          plot.background = element_blank())

    return(gg)
}

## Produces Figure 3 (AUC plots for single antibodies) with IgM subset
##' @param return_data  - return data rather than plot
##' @param subpop - subpopulation to restrict to (shuold be passed as an equation eg., quo(age >= 5))
##' @return dataframe or plot
##' @author Andrew A
make_figure_3 <- function(return_data=FALSE,
                          subpop=NULL){

    ## and let's do another one with only the data where we have IgM LPS
    ## and for IgM
    dat <- load_data_for_paper(long=FALSE,include_igm = TRUE) %>%
    mutate(maxvib=pmax(vibinab,vibogaw))

    if(!is.null(subpop)){
        print("restricting to subpop")
        dat <- dat %>% filter(!!subpop)
    }


    biomarkers_vec <- c('maxvib','pgctxb_s','ptctxa_s','pglsp_s','plpsa_s','plpsm_s','pmctxb_s')
    time_windows <- c('inf_10','inf_45','inf_100','inf_200','inf_365')
    bm_tm <- expand.grid(biomarkers_vec,time_windows)

    single_bm_aucs <- map2(as.character(bm_tm[,1]),
                               as.character(bm_tm[,2]),
                               function(x,y)
                                   cv_titer_repeated_measure_auc(dat,20,x,y) %>% unlist) %>%
    do.call(rbind,.) %>%
    as.data.frame %>%
    mutate(bm=bm_tm[,1],tm=bm_tm[,2]) %>%
      mutate(cvAUC=cvAUC*100,
             ci1=ci1*100,
             ci2=ci2*100)


    ## standarizing colors
    ## should make this into a function at some point
    mycols <- brewer.pal(8,"Dark2")
    names(mycols) <- c('vibinab','vibogaw',
                         'pglsp_s','pgctxb_s',
                         'plpsa_s','ptctxa_s',
                         'plpsm_s','pmctxb_s')
    my_scale <- scale_colour_manual(name = "assay",values = mycols,
                                    labels=c("Vibriocidals",
                                             "Anti-LPS IgG","Anti-CTB IgG",
                                             "Anti-LPS IgA","Anti-CTB IgA",
                                             "Anti-LPS IgM","Anti-CTB IgM"))


    single_bm_aucs %>%
    mutate(bm = as.character(bm),
           bm = if_else(bm == "maxvib","vibogaw",bm),
           bm = keep_assay_color_consistent(bm)) %>%
    ggplot(aes(x=tm,y=cvAUC,color=bm)) +
    geom_pointrange(aes(ymin = ci1, ymax = ci2),position=position_dodge(width = .3)) +
    theme_bw() +
    my_scale +
    theme(legend.position="top") +
    xlab("Infection Time Window") +
    ylab("Cross-Validated Area Under the Curve") +
    scale_x_discrete(labels=c("10-day","45-day","100-day","200-day","365-day"))+
    guides(color=guide_legend(title=NULL)) -> single_marker_auc_gg

    if(return_data) {return(single_bm_aucs)}
    else {return(single_marker_auc_gg)}

}

## makes AUC estimates for single biomarkers
## NOTE THAT THIS FIGURE 3 has been replaced during revisions
##' but is still used for s11
#' @param subpop - used for making estimates for subpoulations but need to provide a quo object (e.g., for filtering for people with age <5 should be quo(age <5)
#' @param return_data - should we return figure or data
make_figure_3_OLD <- function(return_data=FALSE,
                              subpop=NULL){

    dat <- load_data_for_paper(long=FALSE) %>% mutate(maxvib=pmax(vibinab,vibogaw))

    if(!is.null(subpop)){
        print("restricting to subpop")
        dat <- dat %>% filter(!!subpop)
    }

    biomarkers_vec <- c('maxvib','pgctxb_s','ptctxa_s','pglsp_s','plpsa_s')
    time_windows <- c('inf_10','inf_45','inf_100','inf_200','inf_365')
    bm_tm <- expand.grid(biomarkers_vec,time_windows)

    single_bm_aucs <- map2(as.character(bm_tm[,1]),
                           as.character(bm_tm[,2]),
                           function(x,y) cv_titer_repeated_measure_auc(dat,20,x,y) %>% unlist) %>%
    do.call(rbind,.) %>%
    as.data.frame %>%
    mutate(bm=bm_tm[,1],tm=bm_tm[,2]) %>%
      mutate(cvAUC=cvAUC*100,
             ci1=ci1*100,
             ci2=ci2*100)

    mycols <- brewer.pal(8,"Dark2")
    names(mycols) <- c('vibinab','vibogaw',
                         'pglsp_s','pgctxb_s',
                         'plpsa_s','ptctxa_s',
                         'plpsm_s','pmctxb_s')
    my_scale <- scale_colour_manual(name = "assay",values = mycols,
                                    labels=c("Vibriocidals",
                                             "Anti-LPS IgG","Anti-CTB IgG",
                                             "Anti-LPS IgA","Anti-CTB IgA",NA,NA))

    single_bm_aucs %>% mutate(
                           bm = as.character(bm),
                           bm = if_else(bm == "maxvib","vibogaw",bm),
                           bm = keep_assay_color_consistent(bm)) %>%
    ggplot(aes(x=tm,y=cvAUC,color=bm)) +
    geom_pointrange(aes(ymin = ci1, ymax = ci2),position=position_dodge(width = .3)) +
    theme_bw() +
    my_scale +
    theme(legend.position="top") +
    xlab("Infection Time Window") +
    ylab("Cross-Validated Area Under the Curve") +
    scale_x_discrete(labels=c("10-day","45-day","100-day","200-day","365-day"))+
    #scale_color_discrete(labels = c("Vibriocidal", "Anti-CTB IgG","Anti-CTB IgA",
    #                                "Anti-LPS IgG","Anti-LPS IgA"),palette =assay_color_scale() ) +
    guides(color=guide_legend(title=NULL)) -> single_marker_auc_gg

    if(return_data) {return(single_bm_aucs)}
    else {return(single_marker_auc_gg)}

}

## makes figure 4, the figure with ROC curves for RF models
##' @param refit - if we want to refit models
##' @param igm - include IgM subset only?
##' @return
##' @author Andrew A
make_figure_4 <- function(refit=FALSE,igm=FALSE){
    ## get the ROCs and AUC performance for each
    times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')

    if(refit){


        if(igm){
            dat <- load_data_for_paper(long=FALSE,include_igm=TRUE) %>%
            filter(!is.na(plpsm_s) & !is.na(pmctxb_s))

            full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                              'pglsp_s','pgctxb_s','pmctxb_s','plpsm_s',
                              'o_group','sex','wt','age',sep="+")


        } else {
            dat <- load_data_for_paper(long=FALSE)

            full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                              'pglsp_s','pgctxb_s','o_group','sex','wt','age',sep="+")

        }


        ntrees <- 1000
        k <- 20

        set.seed(262340)
        cv_roc_rfs <- lapply(1:length(times),
                             function(t){
                                 my_formula <- paste0("factor(",times[t],") ~",full_set) %>% as.formula
                                 print(my_formula)
                                 cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)
                             })

        print("saving CV models list to generated_data/cv_roc_rfs_full.rds")

        if(igm){
            saveRDS(cv_roc_rfs,"generated_data/cv_roc_rfs_full_igm.rds")
        } else {
            saveRDS(cv_roc_rfs,"generated_data/cv_roc_rfs_full.rds")
        }

    } else {

        if(igm){
            cv_roc_rfs <- readRDS("generated_data/cv_roc_rfs_full_igm.rds")
        } else {
            cv_roc_rfs <- readRDS("generated_data/cv_roc_rfs_full.rds")

        }
    }

    gg <- vector("list",length=5)
    for (i in 1:5){
        preds <- cv_roc_rfs[[i]][[3]]
        truths <- cv_roc_rfs[[i]][[2]]
        imps <- cv_roc_rfs[[i]][[4]]

        days <- str_extract(times[i],"[0-9]+")
        gg[[i]] <- make_roc_varimp_plot(preds,split(truths$truth,truths$fold),imps,
                                        my_title=paste0(days,"-day Infection Window"),
                                        sub=TRUE,ribbon=FALSE,annot_auc=FALSE)
    }

    ## return multiplot object
    multiplot(gg[[1]] + annotate("text",x=0,y=.99,label="A",size=9),
              gg[[3]] + annotate("text",x=0,y=.99,label="C",size=9),
              gg[[2]] + annotate("text",x=0,y=.99,label="B",size=9),
              gg[[5]] + annotate("text",x=0,y=.99,label="D",size=9),
              cols=2)

}

## makes ROC curves for paxvax dataset
##' @author Andrew A
make_figure_5 <- function(){

    set.seed(453453)

    ## load models
    forest200_two <- readRDS("generated_data/forest200reduced.rds")
    forest100_two <- readRDS("generated_data/forest100reduced.rds")
    forest45_two <- readRDS("generated_data/forest45reduced.rds")

    ## load pvdat
    pvdat <- read_csv("data/public_challenge_data.csv") %>%
    filter(!is.na(titer)) %>%
    distinct %>%
    mutate(inf_45=ifelse(day>4&day<=45,1,0),
           inf_100=ifelse(day>4&day<=100,1,0),
           inf_200=ifelse(day>4&day<=200,1,0),
           inf_365=ifelse(day>4&day<=365,1,0),
           titer=ifelse(titer==10 & assay %in% c('vibogaw','vibinab'),5,titer)) %>%
    mutate(dia=cut(dia_vol,c(-Inf,0,3001,Inf),right = TRUE,labels=c("None","Mild","Severe")))


    ## make dataset in a format for the RF models
    rfdat <- pvdat %>% spread(key=assay,value=titer) %>%
    filter(!is.na(vibogaw),!is.na(pgctxb_s))

    ## predicting on all observations
    preds_45 <- predict(forest45_two,newdata = rfdat %>% select(inf_45,vibogaw,pgctxb_s),predict.all=TRUE)
    votes_45 <- apply(preds_45$individual,1,function(x) mean(x==1))
    perf_auc_45 <- prediction(votes_45,rfdat%>% select(inf_45) %>% unlist) %>%
    performance(.,"auc")
    auc_45 <- perf_auc_45@y.values[[1]]
    perf_roc_45 <- prediction(votes_45,rfdat%>% select(inf_45) %>% unlist) %>% performance(.,"tpr","fpr")

    df_45 <- data.frame('fpr'=perf_roc_45@x.values[[1]],'tpr'=perf_roc_45@y.values[[1]],
                        'prob'=perf_roc_45@alpha.values[[1]])

    ## 100
    preds_100 <- predict(forest100_two,newdata = rfdat %>% select(inf_100,vibogaw,pgctxb_s),predict.all=TRUE)
    votes_100 <- apply(preds_100$individual,1,function(x) mean(x==1))
    perf_auc_100 <- prediction(votes_100,rfdat%>% select(inf_100) %>% unlist) %>%
    performance(.,"auc")
    auc_100 <- perf_auc_100@y.values[[1]]
    perf_roc_100 <- prediction(votes_100,rfdat%>% select(inf_100) %>% unlist) %>% performance(.,"tpr","fpr")

    df_100 <- data.frame('fpr'=perf_roc_100@x.values[[1]],'tpr'=perf_roc_100@y.values[[1]],
                         'prob'=perf_roc_100@alpha.values[[1]])

    ## 200
    preds_200 <- predict(forest200_two,newdata = rfdat %>% select(inf_200,vibogaw,pgctxb_s),predict.all=TRUE)
    votes_200 <- apply(preds_200$individual,1,function(x) mean(x==1))
    perf_auc_200 <- prediction(votes_200,rfdat%>% select(inf_200) %>% unlist) %>%
    performance(.,"auc")
    auc_200 <- perf_auc_200@y.values[[1]]
    perf_roc_200 <- prediction(votes_200,rfdat%>% select(inf_200) %>% unlist) %>% performance(.,"tpr","fpr")

    df_200 <- data.frame('fpr'=perf_roc_200@x.values[[1]],'tpr'=perf_roc_200@y.values[[1]],
                         'prob'=perf_roc_200@alpha.values[[1]])


    roc_combined <- bind_rows(data.frame(tw=200,df_200),
                              data.frame(tw=100,df_100),
                              data.frame(tw=45,df_45))

    roc_combined %>%
    mutate(tw=factor(tw,labels=c('45-days','100-days','200-days'))) %>%
    ggplot(aes(x=fpr,y=tpr)) +
    geom_line(aes(color=tw)) +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    annotate("text", x = 0.2, y = 1.0, size=4,alpha=.75,label=paste0("AUC[200]==",sprintf("%.1f",round(auc_200*100,1))),parse=TRUE)+
    annotate("text", x = 0.35, y = .78, size=4,alpha=.75,label = paste0("AUC[100]==",sprintf("%.1f",round(auc_100*100,1))),parse=TRUE)+
    annotate("text", x = 0.22, y = .6, size=4,alpha=.75,label = paste0("AUC[45]==",sprintf("%.1f",round(auc_45*100,1))),parse=TRUE) +
    theme(legend.position = c(.75,.21)) +
    scale_color_brewer(palette="Paired") +
    labs(color="Infection Window") +
    geom_abline(intercept = 0,alpha=.3,lty=2)-> gg

    return(gg)

}

##' makes Figure with results from epidemic
##' illustrative example
##' @title
##' @param rerun - TRUE if we want to re-run all sims (this takes a long time)
##' @return
##' @author Andrew A
make_figure_6 <- function(rerun=FALSE){

    ## exposure to case ratio
    inf_to_case <- c(.5,seq(1,16,by=3))

    ws_epi <- read_csv('data/wau_shilluck_2014.csv') %>%
    mutate(day=(admission_date - as.Date('2014-06-30')) %>% as.numeric)

    N <- 39000

    ## delay to sero-survey
    serosurvey_delay <- 21

    if(rerun){

        my_forest <- readRDS("generated_data/forest365full.rds")

        ## here we will generate the full population profiles
        ## for each realization of the epidemic

        ## realizations per epi curve
        sims_per_epi <- 100
        seros_per_epi <- 100

        ws_epi <- ws_epi %>% mutate(infected=cases,
                                    time_since_inf=max(day)+serosurvey_delay-day) %>%
        data.frame

        ## now for each infection to case ratio and sims_per_epi$cases
        ## generate a population profile and save
        make_infected_times <- function(inf_to_case,ws_epi){
            c(rep(ws_epi$time_since_inf,ws_epi$infected*inf_to_case),
              rep(Inf,N-sum(ws_epi$infected*inf_to_case)))
        }

        ## generates and sames them
        set.seed(5635)
        full_profiles <- purrr::map(inf_to_case,~ generate_epidemic_pop_profiles(make_infected_times(.x,ws_epi),
                                                                                 sims_per_epi,inf_to_case_label=.x))

        saveRDS(full_profiles,"generated_data/wau_shilluck_sim_fullprofs.rds")

        ws_sims_500 <- purrr::map(full_profiles,~purrr::map(.,~ sapply(1:seros_per_epi,function(x) one_sim_singlemodel_v2(500,.,my_forest),simplify=F))) %>%
        flatten %>% flatten %>%
        do.call('rbind',.) %>%
        data.frame

        ws_sims_1000 <- purrr::map(full_profiles,~purrr::map(.,~ sapply(1:seros_per_epi,function(x) one_sim_singlemodel_v2(1000,.,my_forest),simplify=F))) %>%
        flatten %>% flatten %>%
        do.call('rbind',.) %>%
        data.frame

        ws_sims_3000 <- purrr::map(full_profiles,~purrr::map(.,~ sapply(1:seros_per_epi,function(x) one_sim_singlemodel_v2(3000,.,my_forest),simplify=F))) %>%
        flatten %>% flatten %>%
        do.call('rbind',.) %>%
        data.frame

        ws_sims <- bind_rows(ws_sims_500,ws_sims_1000,ws_sims_3000)

        saveRDS(ws_sims,"generated_data/wau_shilluck_sim_overview.rds")

    } else {

        ws_sims <- readRDS("generated_data/wau_shilluck_sim_overview.rds") %>%
        mutate(truth = true_positive_N/39000)
    }

    ## correct vibriocidals
    ## with 200 day sens and spec
    ws_sims <- ws_sims %>% mutate(
                               #truth = true_positive_N/N,
                               cum_inc_vib_cor =
                               pmax(correct_sero_misclass(positive_vib,n-positive_vib,
                                                          sens=0.824,spec=0.81)/n,0),
                               ae_rf = 100*abs(truth - cum_inc_rf)/truth,
                               ae_vib = 100*abs(truth - cum_inc_vib_cor)/truth)

    main_plot <- ws_sims %>% filter(n==500) %>% ggplot() +
    geom_jitter(aes(x=truth*1000,y=cum_inc_rf*1000),
                alpha=.03,color="#3182bd",width=0.03,height=0,size=.5) +
    geom_boxplot(aes(x=truth*1000,y=cum_inc_rf*1000,
                     group=cut_width(truth*1000,5)),
                 alpha=.3,outlier.alpha=.1,lwd=.2,outlier.size=-1) +
    geom_abline(intercept=0,slope=1,lty=2,alpha=.25) +
    scale_x_continuous(trans="log1p",limits=c(10,800),breaks=c(1,5,10,20,50,100,200,400,800)) +
    scale_y_continuous(trans="log1p",limits=c(0,800),breaks=c(1,5,10,20,50,100,200,400,800)) +
    xlab("simulated attack rate (per 1,000)") +
    ylab("estimated attack rate (per 1,000) \n random forest model") +
    theme(axis.title=element_text(size=8)) +
    annotate("text",x=10,y=700,label="B",size=7,fontface =2) +
    annotate("text",x=unique(1000*ws_sims$truth),y=1,label=paste0(inf_to_case,":1"),size=3,angle=90)

    vib_plot <- ws_sims %>%
    filter(n==500) %>%
    ggplot() +
    geom_jitter(aes(x=truth*1000,y=cum_inc_vib_cor*1000),
                alpha=.03,color="#3182bd",width=.03,height=0,size=.5) +
    geom_boxplot(aes(x=truth*1000,y=cum_inc_vib_cor*1000,
                     group=cut_width(truth*1000,5)),
                 alpha=.3,outlier.alpha=.1,lwd=.2,outlier.size = -1) +
    geom_abline(intercept=0,slope=1,lty=2,alpha=.25) +
    scale_x_continuous(trans="log1p",limits=c(10,800),breaks=c(1,5,10,20,50,100,200,400,800)) +
    scale_y_continuous(trans="log1p",limits=c(0,800),breaks=c(1,5,10,20,50,100,200,400,800)) +
    xlab("simulated attack rate (per 1,000)") +
    ylab("estimated attack rate (per 1,000) \n corrected vibriocidal test") +
    theme(axis.title=element_text(size=8)) +
    annotate("text",x=10,y=700,label="C",size=7,fontface =2)

    epi_curves <- data.frame(day=c(ws_epi$day,48),
                             cases=c(ws_epi$cases,0))

    epi_curve <- epi_curves %>%
    ggplot(aes(x=day,y=cases)) +
    geom_bar(
        stat='identity',
        position=position_stack(reverse = TRUE),
        alpha=.9) +
    scale_fill_brewer(palette = "Set2") +
    geom_rect(xmin=66.5,xmax=69.5,ymin=-10,ymax=200,inherit.aes=FALSE,fill="#fc9272") +
    ylab("daily cholera cases") +
    xlab('epidemic day') +
    scale_x_continuous(breaks=c(0,seq(10,70,by=10)),limits=c(0,70)) +
    theme(legend.position = c(0.8, 0.7)) +
    annotate("text",x=0,y=120,label="A",size=7,fontface =2) +
        theme(axis.title=element_text(size=8)) +
    annotate("text",x=67.5,y=50,label="serosurvey",angle=90)

    multiplot(epi_curve,main_plot,vib_plot,cols=2,layout=matrix(c(1,1,2,3),nrow=2,byrow=TRUE))
}

## makes table 1
##' @param target_fn - file name for writing table 1 (use .docx)
##' @return
##' @author Andrew A
make_table_1 <- function(target_fn="table_1.docx"){
    dat <- load_data_for_paper(long=FALSE,
                               extra_vars=c("dehyd")) %>% group_by(id)
    ## sample size
    samp_size <- dat %>% summarize(type=ifelse(cx[1]==9,"contact","case")) %>%
    group_by(type) %>% summarize(n())
    r1 <- data.frame(field="Number of Participants",
                     cases=samp_size[1,2] %>% as.character,
                     contacts=samp_size[2,2] %>% as.character)

    ## age
    age <- dat %>%
    summarize(type=ifelse(cx[1]==9,"contact","case"),age=age[1]) %>%
    group_by(type) %>%
    dplyr::summarize(q50=median(age),q25=round(quantile(age,p=.25),1),q75=round(quantile(age,p=.75),1)) %>%
    mutate(med_age = paste0(q50," (",q25,"-",q75,")")) %>% select(med_age)
    r2 <- data.frame(field="Age of Participants",
                     cases=age[1,1] %>% unname %>% as.character,
                     contacts=age[2,1] %>% unname %>% as.character)

    ## sex
    sex <- dat %>%
    summarize(type=ifelse(cx[1]==9,"contact","case"),sex=sex[1]) %>%
    group_by(type) %>%
    dplyr::summarize(sex=mean(sex))
    r3 <- data.frame(field="Male, %",
                     cases=round(100*sex[1,2],1) %>% unname %>% as.character,
                     contacts=round(100*sex[2,2],1) %>% unname %>% as.character)

    ## O blood group
    obg <- dat %>%
    summarize(type=ifelse(cx[1]==9,"contact","case"),o_group=o_group[1]) %>%
    group_by(type) %>%
    dplyr::summarize(o_group=mean(o_group))
    r4 <- data.frame(field="O blood group, %",
                     cases=round(100*obg[1,2],1) %>% unname %>% as.character,
                     contacts=round(100*(obg[2,2]),1) %>% unname %>% as.character)

    ## Ogawa
    ogawa <- dat %>%
    summarize(type=ifelse(cx[1]==9,"contact","case"),cx=cx[1]) %>%
    group_by(type) %>%
    dplyr::summarize(o_group=mean(cx==1))

    r5 <- data.frame(field="V cholerae O1 Ogawa isolated, %",
                     cases=round(100*ogawa[1,2],1) %>% unname %>% as.character,
                     contacts="-")

    ## dehydration
    dhyd <- dat %>%
    summarize(type=ifelse(cx[1]==9,"contact","case"),dehyd=dehyd[1]) %>%
    group_by(type) %>%
    dplyr::summarize(dehyd=mean(dehyd==3))

    r6 <- data.frame(field="Severely Dehydrated at Admission, %",
                     cases=round(100*dhyd[1,2],1) %>% unname %>% as.character,
                     contacts="-")

    my_tab <- bind_rows(r1,r2,r3,r4,r5,r6) %>% regulartable %>% autofit()

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)


}

## makes table with AUC overview of different RF models
##' @title
##' @param refit
##' @param remove_exposed_contacts
##' @param igm igm subset?
##' @param target_fn filename for output
##' @return
##' @author Andrew A
make_table_2 <- function(refit=FALSE,
                         remove_exposed_contacts = FALSE,
                         igm = FALSE,
                         target_fn = "table_2.docx"){

    if(remove_exposed_contacts & !refit) {
        stop("if you want to do this with remove_exposed_contacts it needs to be refit")
    }

    if(refit){

        if(igm){
            dat <- load_data_for_paper(long=FALSE,include_igm=TRUE) %>%
            filter(!is.na(plpsm_s) & !is.na(pmctxb_s))

            full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                              'pglsp_s','pgctxb_s','pmctxb_s','plpsm_s',
                              'o_group','sex','wt','age',sep="+")


        } else {
            dat <- load_data_for_paper(long=FALSE)

            full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                              'pglsp_s','pgctxb_s','o_group','sex','wt','age',sep="+")

        }

        if(remove_exposed_contacts){
            print("fitting JAGS model for mixture distribution")

            exlc_ids <- id_clusters() %>% filter(cluster==2) %>% select(id) %>% unlist

            dat <- dat %>% filter(!id %in% exlc_ids)

        }

        ntrees <- 1000
        k <- 20
        #full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
        #                  'pglsp_s','pgctxb_s','o_group','sex','wt','age',sep="+")
        times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')

        ## --------------------- ##
        ## First for full models ##
        ## --------------------- ##

        set.seed(262340)
        cv_roc_rfs <- lapply(1:length(times),
                             function(t){
                                 my_formula <- paste0("factor(",times[t],") ~",full_set) %>% as.formula
                                 print(my_formula)
                                 cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)
                             })

        if(!remove_exposed_contacts){

            if(igm){

                print("saving CV models list to generated_data/cv_roc_rfs_igm_full.rds")
                saveRDS(cv_roc_rfs,file="generated_data/cv_roc_rfs_igm_full.rds")
                suffix <- "igm_full"

            } else {

                print("saving CV models list to generated_data/cv_roc_rfs_full.rds")
                saveRDS(cv_roc_rfs,file="generated_data/cv_roc_rfs_full.rds")
                suffix <- "full"

            }


            print("saving individual full models as generated_data/forest_*")
            set.seed(675430)
            for (t in seq_along(times)){
                my_formula <- paste0("factor(",times[t],") ~",full_set) %>% as.formula

                my_forest <- randomForest(my_formula,
                                          ntree=ntrees,
                                          data=dat,
                                          importance=FALSE)
                day <- str_extract(times[t],"[0-9]+")

                saveRDS(my_forest,file=paste0("generated_data/forest",day,suffix,".rds"))
            }
        }

        ## --------------------------------- ##
        ## now for reduced 2 biomarker model ##
        ## --------------------------------- ##

        set.seed(262341)
        cv_roc_rfs_reduced <- lapply(1:length(times),
                                     function(t){

                                         my_set <- paste('vibogaw','pgctxb_s',sep="+")

                                         if(t==1){
                                             my_set <- paste('vibogaw','ptctxa_s',sep="+")
                                         }

                                         my_formula <- paste0("factor(",times[t],") ~",my_set) %>% as.formula
                                         print(my_formula)

                                         cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)#,
                                         ##cutoff=c(reduced_cutoffs[t],1-reduced_cutoffs[t]))
                                     })

        if(!remove_exposed_contacts){

            if(igm){
                print("saving CV models list to generated_data/cv_roc_rfs_reduced.rds")
                saveRDS(cv_roc_rfs_reduced,file="generated_data/cv_roc_rfs_reduced.rds")
                suffix <- "igm_reduced.rds"
            } else {

                print("saving CV models list to generated_data/cv_roc_rfs_igm_reduced.rds")
                saveRDS(cv_roc_rfs_reduced,file="generated_data/cv_roc_rfs_igm_reduced.rds")
                suffix <- "reduced.rds"
            }

            ## now save these models
            print("saving 2 biomarker models")
            ## save reduced MODELS
            set.seed(67430)
            for (t in seq_along(times)){
                my_set <- paste('vibogaw','pgctxb_s',sep="+")
                if(t==1){
                    my_set <- paste('vibogaw','ptctxa_s',sep="+")
                }
                my_formula <- paste0("factor(",times[t],") ~",my_set) %>% as.formula

                my_forest <- randomForest(my_formula,
                                          ntree=ntrees,
                                          data=dat,
                                          importance=FALSE)#,
                                        #cutoff=c(reduced_cutoffs[t],1-reduced_cutoffs[t]))

                day <- str_extract(times[t],"[0-9]+")

                saveRDS(my_forest,file=paste0("generated_data/forest",day,suffix))
            }

        }

        ## ------------------------- ##
        ## Now for ELISA only models ##
        ## ------------------------- ##
        set.seed(264341)

        if(igm){

            ELISA_set <- paste('plpsa_s','ptctxa_s',
                               'pglsp_s','pgctxb_s',
                               'plpsm_s','pmctxb_s',sep="+")

        } else {

            ELISA_set <- paste('plpsa_s','ptctxa_s',
                               'pglsp_s','pgctxb_s',sep="+")

        }

        cv_roc_rfs_ELISA <- lapply(1:length(times),
                                   function(t){
                                       my_formula <- paste0("factor(",times[t],") ~",ELISA_set) %>% as.formula
                                       print(my_formula)
                                       cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)
                                   })

        if(!remove_exposed_contacts){

            if(igm){
                print("saving CV models list to generated_data/cv_roc_rfs_igm_ELISA.rds")
                saveRDS(cv_roc_rfs_ELISA,file="generated_data/cv_roc_rfs_igm_ELISA.rds")
                suffix <- "igm_ELISA.rds"
            } else {
                print("saving CV models list to generated_data/cv_roc_rfs_ELISA.rds")
                saveRDS(cv_roc_rfs_ELISA,file="generated_data/cv_roc_rfs_ELISA.rds")
                suffix <- "ELISA.rds"
            }

            print("saving ELISA models")

            set.seed(67430)
            for (t in seq_along(times)){
                my_formula <- paste0("factor(",times[t],") ~",ELISA_set) %>% as.formula

                my_forest <- randomForest(my_formula,
                                          ntree=ntrees,
                                          data=dat,
                                          importance=FALSE)

                day <- str_extract(times[t],"[0-9]+")
                saveRDS(my_forest,file=paste0("generated_data/forest",day,suffix))
            }
        }

    } else {

        if(igm){
            cv_roc_rfs <- readRDS("generated_data/cv_roc_rfs_igm_full.rds")
            cv_roc_rfs_reduced <- readRDS("generated_data/cv_roc_rfs_igm_reduced.rds")
            cv_roc_rfs_ELISA <- readRDS("generated_data/cv_roc_rfs_igm_ELISA.rds")
        } else {
            cv_roc_rfs <- readRDS("generated_data/cv_roc_rfs_full.rds")
            cv_roc_rfs_reduced <- readRDS("generated_data/cv_roc_rfs_reduced.rds")
            cv_roc_rfs_ELISA <- readRDS("generated_data/cv_roc_rfs_ELISA.rds")

        }

    }

    ## have a look at CV AUCs
    overview_full <- map_dfr(cv_roc_rfs,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_full) <- "full"


    ## have a look at performance of reduced model
    overview_reduced <- map_dfr(cv_roc_rfs_reduced,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_reduced) <- "reduced"

    overview_ELISA <- map_dfr(cv_roc_rfs_ELISA,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_ELISA) <- "ELISA"

    overview_tab <- (bind_rows(overview_full,overview_reduced,overview_ELISA)) %>% set_rownames(c("full model","two biomarker","ELISA only"))

    my_tab <- regulartable(bind_cols(data.frame("Model"=rownames(overview_tab)),overview_tab)) %>%
    add_header(Model = "",
               `10-days`="Infection Time Window",
               `45-days`="Infection Time Window",
               `100-days`="Infection Time Window",
               `200-days`="Infection Time Window",
               `365-days`="Infection Time Window"
               ) %>% merge_h(part = "header") %>%
    align(align="center",part="all") %>%
    autofit

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

    print(my_tab)

}


##' makes table S1
##' AUCs for single biomarkers and 95% CIs (based on influence curves)
##' @param target_fn filename where table should be saved (should be .docx)
##' @return
##' @author Andrew A
make_table_s1 <- function(target_fn = "table_s1.docx"){
    my_dat <- make_figure_3(return_data = TRUE)
    my_dat <- my_dat %>%
    mutate(
        "cvAUC (95% CI)" = paste0(round(cvAUC,1),"% \n(",round(ci1,1),"-",round(ci2,1),")"),
        "Time Window" = clean_inf_wind_labels(tm),
        "Biomarker" = clean_biomarker_labels(bm)) %>%
    select(`cvAUC (95% CI)`,`Time Window`, `Biomarker`) %>% spread(`Biomarker`,`cvAUC (95% CI)`)

    my_tab <- regulartable(my_dat) %>% width(width=1.4) %>% align(align="center",part="all") %>%
    bold(j = 1, bold = TRUE, part = "body")

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

    print(my_tab)
}


##' Makes Table S2, which contains estimates of sens and spec
##' for single biomarkers
##' @title
##' @param target_fn
##' @param rerun
##' @return
##' @author Andrew A
make_table_s2 <- function(target_fn = "table_s2.docx",
                          rerun = TRUE){
    if(rerun){
        set.seed(6548)
        dat <- load_data_for_paper(long=FALSE,include_igm=TRUE)
        dat <- dat %>%  mutate(maxvib=pmax(vibinab,vibogaw))
        times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')

        ## go through each biomarker
        print("vibriocidals")
        vibs_perf <- get_single_bm_sens_spec(dat,'maxvib')
        vib_out <- vibs_perf[[1]] %>% mutate(bm = "Vibriocidal")

        print("CTB IgG")
        ctb_igg_perf <- get_single_bm_sens_spec(dat,'pgctxb_s')
        ctb_igg_out <- ctb_igg_perf[[1]] %>% mutate(bm = "CTB IgG")

        print("CTB IgA")
        ctb_iga_perf <- get_single_bm_sens_spec(dat,'ptctxa_s')
        ctb_iga_out <- ctb_iga_perf[[1]] %>% mutate(bm = "CTB IgA")

        print("CTB IgM")
        ctb_igm_perf <- get_single_bm_sens_spec(dat,'pmctxb_s')
        ctb_igm_out <- ctb_igm_perf[[1]] %>% mutate(bm = "CTB IgM")

        print("LPS IgG")
        lps_igg_perf <- get_single_bm_sens_spec(dat,'pglsp_s')
        lps_igg_out <- lps_igg_perf[[1]] %>% mutate(bm = "LPS IgG")

        print("LPS IgA")
        lps_iga_perf <- get_single_bm_sens_spec(dat,'plpsa_s')
        lps_iga_out <- lps_iga_perf[[1]] %>% mutate(bm = "LPS IgA")

        print("LPS IgM")
        lps_igm_perf <- get_single_bm_sens_spec(dat,'plpsm_s')
        lps_igm_out <- lps_igm_perf[[1]] %>% mutate(bm = "LPS IgM")


        rc <- bind_rows(vib_out,ctb_igg_out,ctb_igm_out,ctb_iga_out,lps_igg_out,lps_igm_out,lps_iga_out)
        saveRDS(rc,"generated_data/single_biomarker_performance.rds")

    } else {
        ## reading saved performance data
        rc <- readRDS("generated_data/single_biomarker_performance.rds")

    }

    my_tab <- rc %>% ungroup %>%
    filter(variable %in% c('sens_test','spec_test')) %>%
    mutate(time=recode(time,
                       inf_10='10 days',
                       inf_45='45 days',
                       inf_100='100 days',
                       inf_200='200 days',
                       inf_365='365 days'
                       )) %>%
    group_by(bm,time) %>%
    dplyr::summarize(sensitivity=sprintf("%.1f%% (%.1f-%.1f)",mean[1]*100,p2.5[1]*100,p97.5[1]*100),
                     specificity=sprintf("%.1f%% (%.1f-%.1f)",mean[2]*100,p2.5[2]*100,p97.5[2]*100),
                     `modal titer`=round(modal_titer[1],1))%>%
    rename('Biomarker' = bm,
           'Time Window' = time) %>%
    regulartable() %>%
    merge_v(j="Biomarker") %>%
    width(width=1.4) %>%
    align(align="center",part="all") %>%
    bold(j = 1, bold = TRUE, part = "body")


    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

    print(my_tab)
}

##' makes table S3
##' table S2 but with IgM subset
##' @return html table (and writes to word doc)
##' @author Andrew A
make_table_s3 <- function(){
    make_table_2(refit=TRUE, igm=TRUE,target_fn = "table_s3.docx")
}

##' Makes table S4
##' performance table for O-negative particpants only
##' @return
##' @author Andrew A
make_table_s4 <- function(){
    make_aucs_oneg_only(target_fn = "table_s4.docx")
}

##' Makes table S5
##' Comparison of AUC models classifying people in adjacent
##' time windows
## NOTE: doing 10-fold cross valdation, not 20 given small sample size
##' @param target_fn
##' @return
##' @author Andrew A
make_table_s5 <- function(target_fn="table_s5.docx"){

    ntrees <- 1000
    full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                      'pglsp_s','pgctxb_s','o_group','sex','wt','age',sep="+")
    times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')

    dat_temp_sig <- dat %>% mutate(inf_10_45 = case_when(
                                       inf_10 ==1 ~ 1,
                                       (inf_45 == 1 & inf_10 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_10_100 = case_when(
                                       inf_10 ==1 ~ 1,
                                       (inf_100 == 1 & inf_45 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_10_200 = case_when(
                                       inf_10 ==1 ~ 1,
                                       (inf_200 == 1 & inf_100 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_10_365 = case_when(
                                       inf_10 ==1 ~ 1,
                                       (inf_365 == 1 & inf_200 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_45_100 = case_when(
                                       inf_45 ==1 & inf_10 == 0 ~ 1,
                                       (inf_100 == 1 & inf_45 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_45_200 = case_when(
                                       inf_45 == 1 & inf_10 == 0 ~ 1,
                                       (inf_200 == 1 & inf_100 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_45_365 = case_when(
                                       inf_45 ==1 & inf_10 == 0 ~ 1,
                                       (inf_365 == 1 & inf_200 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_100_200 = case_when(
                                       inf_100 == 1 & inf_45 == 0 ~ 1,
                                       (inf_200 == 1 & inf_100 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_100_365 = case_when(
                                       inf_100 ==1 & inf_45 == 0 ~ 1,
                                       (inf_365 == 1 & inf_100 == 0) ~ 0,
                                       TRUE ~ NA_real_),
                                   inf_200_365 = case_when(
                                       inf_200 ==1 & inf_100 == 0 ~ 1,
                                       (inf_365 == 1 & inf_200 == 0) ~ 0,
                                       TRUE ~ NA_real_)
                                   )

    set.seed(356093)
    k <- 10
    sub_times_10 <-c(45,100,200,365)
    cv_roc_rfs_10 <- lapply(1:length(sub_times_10),
                            function(t){
                                tmp_time <- paste0("inf_10_",sub_times_10[t])
                                my_formula <- paste0("factor(",tmp_time,") ~",full_set) %>% as.formula
                                print(my_formula)
                                cv_rf_model(dat_temp_sig %>% filter(!is.na(UQ(rlang::sym(tmp_time)))),
                                            my_formula,
                                            k=k,
                                            ntrees=ntrees)
                            })


    sub_times_45 <-c(100,200,365)
    cv_roc_rfs_45 <- lapply(1:length(sub_times_45),
                            function(t){
                                tmp_time <- paste0("inf_45_",sub_times_45[t])
                                my_formula <- paste0("factor(",tmp_time,") ~",full_set) %>% as.formula
                                print(my_formula)
                                cv_rf_model(dat_temp_sig %>% filter(!is.na(UQ(rlang::sym(tmp_time)))),
                                            my_formula,
                                            k=k,
                                            ntrees=ntrees)
                            })

    sub_times_100 <-c(200,365)
    cv_roc_rfs_100 <- lapply(1:length(sub_times_100),
                             function(t){
                                 tmp_time <- paste0("inf_100_",sub_times_100[t])
                                 my_formula <- paste0("factor(",tmp_time,") ~",full_set) %>% as.formula
                                 print(my_formula)
                                 cv_rf_model(dat_temp_sig %>% filter(!is.na(UQ(rlang::sym(tmp_time)))),
                                             my_formula,
                                             k=k,
                                             ntrees=ntrees)
                             })


    sub_times_200 <-c(365)
    k <- 10
    cv_roc_rfs_200 <- lapply(1:length(sub_times_200),
                             function(t){
                                 tmp_time <- paste0("inf_200_",sub_times_200[t])
                                 my_formula <- paste0("factor(",tmp_time,") ~",full_set) %>% as.formula
                                 print(my_formula)
                                 cv_rf_model(dat_temp_sig %>% filter(!is.na(UQ(rlang::sym(tmp_time)))),
                                             my_formula,
                                             k=k,
                                             ntrees=ntrees)
                             })


    ## helper function to make pretty rows for table
    make_pretty_row <- function(cvs,
                                col_names=paste0(c("10-45","45-100","100-200","200-365")," days"),
                                row_name="m10"){
        rc <- map_dfr(cvs,1) %>%
        mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
        select(lab) %>%
        t %>% unname %>% data.frame %>% set_names(col_names)
        rownames(rc) <- row_name
        return(rc)
    }

    overall_temporal_aucs <- rbind.fill(make_pretty_row(cv_roc_rfs_10),
                                        make_pretty_row(cv_roc_rfs_45,
                                                        col_names = paste0(c("45-100","100-200","200-365")," days"),
                                                        row_name = "m45"),
                                        make_pretty_row(cv_roc_rfs_100,
                                                        col_names = paste0(c("100-200","200-365")," days"),
                                                        row_name = "m100"),
                                        make_pretty_row(cv_roc_rfs_200,
                                                        col_names = paste0(c("200-365")," days"),
                                                        row_name = "m200")
                                        )

    overall_temporal_aucs$cases <- c("5-10 days","10-45 days","45-100 days","100-200 days")


    my_tab <- overall_temporal_aucs %>%
    select(cases, `10-45 days`:`200-365 days`) %>%
    regulartable() %>%
    add_header(`10-45 days`="Comparison Group Infection Window",
               `45-100 days`="Comparison Group Infection Window",
               `100-200 days`="Comparison Group Infection Window",
               `200-365 days`="Comparison Group Infection Window",
               cases = "") %>%
    merge_h(part = "header") %>%
    align(align="center",part="header") %>%
    autofit()


    print(my_tab)

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

}

##' makes table S6
##' table of paxvax sens and spec for single biomarkers
##' @param target_fn
##' @return
##' @author Andrew A
make_table_s6 <- function(target_fn="table_s6.docx"){

    pvdat <- read_csv("data/public_challenge_data") %>%
    filter(!is.na(titer))%>%
    distinct %>%
    mutate(inf_45=ifelse(day>4&day<=45,1,0),
           inf_100=ifelse(day>4&day<=100,1,0),
           inf_200=ifelse(day>4&day<=200,1,0),
           inf_365=ifelse(day>4&day<=365,1,0),
           titer=ifelse(titer==10 & assay %in% c('vibogaw','vibinab'),5,titer)) %>%
    mutate(dia=cut(dia_vol,c(-Inf,0,3001,Inf),right = TRUE,labels=c("None","Mild","Severe")))

    ## not including inf_10 here since we dont' have any samples >4 and <10
    vibriodical_thresholds <- c(640,640,320)
    ctx_igg_thresholds <- c(37,37,37)
    ctx_iga_thresholds <- c(13,9,8)

    thresholds <- c(vibriodical_thresholds,
                    ctx_igg_thresholds,
                    ctx_iga_thresholds)

    assay = c('max_vib','pgctxb_s','ptctxa_s')
    inf_tw = c('inf_45','inf_100','inf_200')

    ## look ups
    lu_tab <- expand.grid(assay=assay,tw=inf_tw) %>%
    arrange(assay) %>%
    mutate(thresholds=thresholds)


    estim_sens_spec <- function(dat,assay_name,thresholds){

        x <- filter(dat,assay==assay_name)

        bind_rows(
            data.frame(time=45,sens_spec(x$inf_45,factor(x$titer > thresholds[1],levels=c(FALSE,TRUE)))),
            data.frame(time=100,sens_spec(x$inf_100,factor(x$titer > thresholds[2],levels=c(FALSE,TRUE)))),
            data.frame(time=200,sens_spec(x$inf_200,factor(x$titer > thresholds[3],levels=c(FALSE,TRUE)))))
    }


    pv_sens_spec <- suppressWarnings(
        bind_rows(
            data.frame(Biomarker='Vibriocidal',estim_sens_spec(pvdat,'max_vib',vibriodical_thresholds)),
            data.frame(Biomarker='Anti-CTB IgG',estim_sens_spec(pvdat,'pgctxb_s',ctx_igg_thresholds)),
            data.frame(Biomarker='Anti-CTB IgA',estim_sens_spec(pvdat,'ptctxa_s',ctx_iga_thresholds)))
    )

    pv_sens_spec <- pv_sens_spec %>%
    rename(Sensitivity=sens,Specificity=spec,`Time Window` = time) %>%
    mutate(Sensitivity = round(100* Sensitivity ,1) %>% as.character,
           Specificity=  round(100*Specificity,1) %>% as.character,
           `Time Window` = `Time Window` %>% as.character)

    my_tab <- pv_sens_spec %>% regulartable %>% merge_v(j=1) %>%
    align(align="center",part="all") %>% autofit
    print(my_tab)

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

    return(pv_sens_spec)
}


##' makes Table S7
##' overview of endemic serosims
##' @param rerun rerun or use saved outputs (this takes a long time!)
##' @param target_fn
##' @return
##' @author Andrew A
make_table_s7 <- function(rerun=FALSE,target_fn="table_s7.docx"){

    N <- 25000 ## population size

    if(rerun){

        trajs <-  make_all_trajs(serotypes=c(1,2)) ## only
        data <- load_data_for_paper(long=FALSE)

        ns <- c(500,1000,3000) ## sample size
        endemic_yrs_per_inc <- 300
        seros_per_epidemic <- 100
        incs <- c(2/1000,5/1000,10/1000,50/1000,100/1000)
        ## annual infection incidence

        my_forest <- readRDS(paste0("generated_data/forest",365,"full.rds"))
        time_since_infection <- vector(mode='list', length=length(incs))
        one_over <- vector(mode='list', length=length(incs))
        my_vib_threshold <- 320

        set.seed(63424)

        ## generate endemic_yrs of populations
        ## for all the incidence
        ## list of length(inc) with a data.frame in each
        ## generate pop profiles to rds files
        pop_profile_file_names <- purrr::map(incs,~generate_endemic_pop_profiles(.x,N,endemic_yrs_per_inc) %>% bind_rows) %>%
        bind_rows()

        pop_profile_file_names <- paste0("generated_data/endemic_pop_prof_inc_",rep(incs * 1000,each=endemic_yrs_per_inc),"_",rep(1:endemic_yrs_per_inc,length(incs)),".rds")


        ## now for each incidence level, load pops one by one and generate a file with all simulated
        ## serosurveys in a single dataframe with a new field called sero_survey_num
        sims_overview <- NULL
        for(i in seq_along(incs)){
            inc_data <- lapply(grep(paste0("inc_",incs[i]*1000),pop_profile_file_names,value=TRUE),readRDS)
            print(incs[i]*1000)
            n1 = purrr::map(inc_data, function(x)
                purrr::map(1:seros_per_epidemic, function(y) one_sim_singlemodel_v2(ns[1],x,my_forest)) %>%
                bind_rows) %>%
            bind_rows
            print("finished n1")
            n2 = purrr::map(inc_data, function(x) purrr::map(1:seros_per_epidemic, function(y) one_sim_singlemodel_v2(ns[2],x,my_forest)) %>% bind_rows) %>% bind_rows
            print("finished n2")
            n3 = purrr::map(inc_data, function(x) purrr::map(1:seros_per_epidemic, function(y) one_sim_singlemodel_v2(ns[3],x,my_forest)) %>% bind_rows) %>% bind_rows
            print("finished n3")
            sims_overview <- bind_rows(sims_overview,n1,n2,n3)
        }

        saveRDS(sims_overview,file="generated_data/endemic_sims_overview.rds")

    } else {

        sims_overview <- readRDS("generated_data/endemic_sims_overview.rds")

    }

    ## get corrected incidence
    ## using 365 day numbers from table s2
    out <- sims_overview %>% mutate(
                                 inc = true_positive_N/N,
                                 positive_vib_cor = correct_sero_misclass(
                                     num_pos = positive_vib,
                                     num_neg = n - positive_vib,
                                     sens=.81,
                                     spec=.83),
                                 cum_inc_vib_cor=positive_vib_cor/n)

    my_tab <- out %>%
    group_by(n,inc) %>%
    summarize(
        `Mean Incidence` = (1000*mean(cum_inc_rf)) %>% round(.,1),
        `Mean Absolute Error (%)`=(100*mean(abs(cum_inc_rf-inc))/inc[1]) %>% round(.,1),
        `Bias (%)`=(100*mean(cum_inc_rf-inc)/inc[1]) %>% round(.,1),
        `Mean Incidence ` = (1000*mean(cum_inc_vib_cor)) %>% round(.,1),
        `Mean Absolute Error (%) `=(100*mean(abs(cum_inc_vib_cor-inc))/inc[1]) %>% round(.,1),
        `Bias (%) `=(100*mean(cum_inc_vib_cor-inc)/inc[1]) %>% round(.,1)) %>%
    mutate(inc = round(inc*1000,1)) %>%
    rename(`Sample Size`=n,`Annual Incidence \n per 1,000`=inc) %>%
    regulartable() %>%  merge_v(j=1) %>% align(align="center",part="all") %>%
    set_formatter_type(fmt_double="%.01f") %>%
    add_header(`Mean Incidence` = "Random Forest",
               `Mean Absolute Error (%)` = "Random Forest",
               `Bias (%)` = "Random Forest",
               `Mean Incidence ` = "Vibriocidal, Corrected",
               `Mean Absolute Error (%) ` = "Vibriocidal, Corrected",
               `Bias (%) ` = "Vibriocidal, Corrected" ,
               top = TRUE ) %>%
    merge_h(part = "header") %>%
    autofit()

    print(my_tab)

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

}

##' Makes table S8
##' Table 2 but only with 'unexposed' contacts
##' @return
##' @author Andrew A
make_table_s8 <- function(){
    make_table_2(refit=TRUE,remove_exposed_contacts = TRUE,target_fn = "table_s8.docx")
}


##' makes figure S1 -  enrollement over time plot
##' @return grid.arrange() object
##' @author Andrew A
make_figure_s1 <- function(){

    stop("this function does not work since we did not share the date of enrollement in the study due for privacy reasons")

    dat_f <- load_data_for_paper(full=TRUE,long=FALSE) %>% mutate(datenrol=as.Date(datenrol))


    dat_f %>% filter(is_case==1) %>%
    ggplot(aes(x=lbday2)) +
    geom_histogram(binwidth=1) +
    xlab('days since case presentation') +
    ylab('number of individuals with blood collected') +
    theme(legend.position=c(.9,.85)) +
    scale_y_sqrt(breaks=c(10,50,100,200,300,400)) +
    annotate(geom="text",x=c(2,7,30,90,180,270,365,540,720,900),y=c(300,250,200,200,200,200,200,200,200,200),
             label=paste0("Day ",c(2,7,30,90,180,270,365,540,720,900))) +
    guides(fill=guide_legend(title="")) + ggtitle("Confirmed Cholera Cases") -> gg


    dat_f %>% filter(is_case==0) %>%
    ggplot(aes(x=lbday2)) +
    geom_histogram(binwidth=1) +
    xlab('days since (household) case presentation') +
    ylab('number of individuals with blood collected') +
    theme(legend.position=c(.9,.85)) +
    scale_y_sqrt() +
    annotate(geom="text",x=c(2,7,30),y=c(75,75,75),
             label=paste0("Day ",c(2,7,30))) +
    guides(fill=guide_legend(title="")) + ggtitle("Household Contacts") -> gg2

    return(
        grid.arrange(gg,gg2,layout_matrix = rbind(c(1,1,1),c(2,2,NA)))
    )

}


##' Makes Figure S2
##' Density plots for CTB antibodies with Igms
##' @return gg object
##' @author Andrew A
make_figure_s2 <- function(){
    gg <- make_titer_dist_density_plots(assays=c('pgctxb_s','pmctxb_s','ptctxa_s'))
    gg <- gg + scale_x_continuous(trans="log2") + xlab("log2(titer)") +
        theme(axis.ticks = element_blank(), legend.background = element_blank(),
          legend.key = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), strip.background = element_blank(),
          plot.background = element_blank())

    return(gg)
}


##' Makes Figure S3
##' denisty plot for LPS including IgMs
##' @return gg object
##' @author Andrew A
make_figure_s3 <- function(){
    gg <- make_titer_dist_density_plots(assays=c('pglsp_s','plpsm_s','plpsa_s'))
    gg <- gg + scale_x_continuous(trans="log2") + xlab("log2(titer)") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(),
          legend.key = element_blank(), panel.background = element_blank(),
          panel.border = element_blank(), strip.background = element_blank(),
          plot.background = element_blank())

    return(gg)

}

##' Makes Figure S4
##' scam plot for vibogaw
##' @return
##' @author Andrew A
make_figure_s4 <- function(){
    rc <- scam_and_plot(assay='vibogaw')
    return(rc)
}

##' Makes Figure S5
##' makes scam plot for vibinab
##' @title
##' @return
##' @author Andrew A
make_figure_s5 <- function(){
    rc <- scam_and_plot('vibinab',ymax=12,ybreaks=2:12,ylabs=2^(2:12))
    return(rc)
}

##' Makes Figure S6
##' SCAM plot for IgG CTB
##' @return
##' @author Andrew A
make_figure_s6 <- function(){
    rc <- scam_and_plot('pgctxb_s',
                        ymax=log2(200),
                        ybreaks=1:9,
                        ylabs=2^(1:9))
    return(rc)
}


##' Makes Figure S7
##' SCAM plot for IgA CTB
##' @return
##' @author Andrew A
make_figure_s7 <- function(){
    rc <- scam_and_plot('ptctxa_s',
                        ymax=log2(200),
                        ybreaks=1:9,
                        ylabs=2^(1:9)
                        )
    return(rc)
}

##' Makes Figure S8
##' SCAM plot for IgG LPS
##' @return
##' @author Andrew A
make_figure_s8 <- function(){
    rc <- scam_and_plot('pglsp_s',
                        ymax=log2(200),
                        ymin = log2(8),
                        ybreaks=3:9,
                        ylabs=2^(3:9))

    return(rc)
}

##' Makes Figure S9
##' SCAM plot for IgA LPS
##' @title
##' @return
##' @author Andrew A
make_figure_s9 <- function(){
    rc <- scam_and_plot('plpsa_s',
                        ymax=log2(200),
                        ybreaks=2:9,
                        ylabs=2^(2:9))

    return(rc)
}

##' Makes Figure S10
##' IgM LPS SCAM plot
##' @return list with gg objects and model comparison stats
##' @author Andrew A
make_figure_s10 <- function(){
    rc <- scam_and_plot('plpsm_s',
                        my_igm = TRUE,
                        ymax=8,
                        ymin = 3,
                        ybreaks=3:10,
                        ylabs=2^(3:10))

    return(rc)
}



##' Figure S11
##' makes variant of Figure 3 with panels for u5,o5 and full model
##' uses older version of make_figure_3
##' @return
##' @author Andrew A
make_figure_s11 <- function(){

    o5 <- make_figure_3_OLD(subpop = quo(age >= 5)) + ggtitle("Five years and older")
    u5 <- make_figure_3_OLD(subpop = quo(age < 5)) + ggtitle("Under five years old")
    full <- make_figure_3_OLD() + ggtitle("All ages")

    grid.arrange(full,u5,o5,nrow=3)
}

##' Makes Figure S12
##' Same as Figure 4 (ROC curves for RF model) but
##' with IgM subset
##' @return
##' @author Andrew A
make_figure_s12 <- function(){
    make_figure_4(refit=TRUE,igm=TRUE)
}


##' Makes Figure S13
##' vibriocidal comparison between pv and dhaka
##' @return gg object
##' @author Andrew A
make_figure_s13 <- function(){

    ## load paxvax data
    pvdat <- read_csv("data/public_challenge_data.csv") %>%
    filter(!is.na(titer)) %>%
    distinct %>%
    mutate(inf_45=ifelse(day>4&day<=45,1,0),
           inf_100=ifelse(day>4&day<=100,1,0),
           inf_200=ifelse(day>4&day<=200,1,0),
           inf_365=ifelse(day>4&day<=365,1,0),
           titer=ifelse(titer==10 & assay %in% c('vibogaw','vibinab'),5,titer)) %>%
    mutate(dia=cut(dia_vol,c(-Inf,0,3001,Inf),right = TRUE,labels=c("None","Mild","Severe")))

    ## load Dhaka data
    dat_l_f <- load_data_for_paper(long=TRUE,full=TRUE)
    augmented_dat_l_f <-
    rbind(subset(dat_l_f,(cx %in% c(1,2) & lbday>=2)),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=2),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=7),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=30),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=90),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=180),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=360)
          )

    controls_only <- augmented_dat_l_f %>% filter(cx==9) %>% ungroup %>%  select(assay,value) %>% mutate(population="Dhaka")

    pvdat_simp <- pvdat %>% filter(day==0) %>% mutate(value=titer,population="USA") %>% select(assay,value,population)

    baseline_comp <- bind_rows(controls_only,pvdat_simp)

    gg <- baseline_comp %>%
    filter(assay %in% c('vibinab','vibogaw')) %>%
    ggplot(aes(x=value)) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") +
    geom_density(aes(fill=population),col=0,alpha=.8) +#,adjust=1.0,bw='nrd')+
    labs(colour = "") +
    ylab("Density") + xlab("vibriocidal titer")+
    facet_grid(assay~.,labeller = as_labeller(assay_names),scale='free_y') +
    geom_jitter(aes(color=population,y=.33),height=.01,width=.1,alpha=.1) +
    scale_x_continuous(trans="log2") +
    guides(color=FALSE) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "#ffffffaa", colour = NA))
    return(gg)
}

##' Makes figure S14
##' comparison of Dhaka and NA volunteers
##' @return gg object
##' @author Andrew A
make_figure_s14 <- function(){
    ## load paxvax data
    pvdat <- read_csv("data/public_challenge_data.csv") %>%
    filter(!is.na(titer)) %>%
    distinct %>%
    mutate(inf_45=ifelse(day>4&day<=45,1,0),
           inf_100=ifelse(day>4&day<=100,1,0),
           inf_200=ifelse(day>4&day<=200,1,0),
           inf_365=ifelse(day>4&day<=365,1,0),
           titer=ifelse(titer==10 & assay %in% c('vibogaw','vibinab'),5,titer)) %>%
    mutate(dia=cut(dia_vol,c(-Inf,0,3001,Inf),right = TRUE,labels=c("None","Mild","Severe")))

    ## load Dhaka data
    dat_l_f <- load_data_for_paper(long=TRUE,full=TRUE)
    augmented_dat_l_f <-
    rbind(subset(dat_l_f,(cx %in% c(1,2) & lbday>=2)),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=2),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=7),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=30),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=90),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=180),
          subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=360)
          )

    controls_only <- augmented_dat_l_f %>% filter(cx==9) %>% ungroup %>%  select(assay,value) %>% mutate(population="Dhaka")

    pvdat_simp <- pvdat %>% filter(day==0) %>% mutate(value=titer,population="USA") %>% select(assay,value,population)

    baseline_comp <- bind_rows(controls_only,pvdat_simp)

    gg <- baseline_comp %>%
    filter(assay %in% c('pgctxb_s','ptctxa_s','pmctxb_s')) %>%
    ggplot(aes(x=value)) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") +
    geom_density(aes(fill=population),col=0,alpha=.8) +#,adjust=1.0,bw='nrd')+
    labs(colour = "") +
    ylab("Density") + xlab("titer")+
    facet_grid(assay~.,labeller = as_labeller(assay_names),scale='free_y') +
    geom_jitter(aes(color=population,y=.65),height=.01,width=.1,alpha=.1) +
    scale_x_continuous(trans="log2") +
    guides(color=FALSE) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "#ffffffaa", colour = NA)
          )
    return(gg)
}


##' Makes Figure S15
##' CTB repsonses in NA volunteers
##' @return gg object
##' @author Andrew A
make_figure_s15 <- function(){
    make_ctb_paxvac_plot()
}


##' utility function to cleans labels
##' @param col column passed
##' @return
##' @author Andrew A
clean_inf_wind_labels <- function(col){
    recode(col,
           inf_10 = "10-day",
           inf_45 = "45-day",
           inf_100 = "100-day",
           inf_200 = "200-day",
           inf_365 = "365-day")
    }

##' utility function for cleaning
##' biomarker labels
##' used for tables
##' @param col
##' @return
##' @author Andrew A
clean_biomarker_labels <- function(col){
    recode(col,
           `plpsa_s` = "Anti-LPS IgA",
           `ptctxa_s` = "Anti-CTB IgA",
           `vibinab` = "Inaba Vibriocidal",
           `vibogaw` = "Ogawa Vibriocidal",
           `maxvib` = "Vibriocidal",
           `max_vib` = "Vibriocidal",
           `pglsp_s` = "Anti-LPS IgG",
           `pgctxb_s` = "Anti-CTB IgG",
           `plpsm_s` = "Anti-LPS IgM",
           `pmctxb_s` = "Anti-CTB IgM")
}


##' makes titer density plots like Figure 2 (and sup figs) in manuscript
##' @param assays - whivh assays to include
##' @param exc_lbdays - which lbdays to exclude for plot (note that 0 is not real)
##' @return
##' @author Andrew A
make_titer_dist_density_plots <- function(assays=c('vibinab','vibogaw'),
                                          exc_lbdays = c(0,270,720)){

    inc_igm <- any(c("plpsm_s","pmctxb_s") %in% assays)

    dat_l_f <- load_data_for_paper(long=TRUE,full=TRUE,include_igm = inc_igm)

    augmented_dat_l_f <-
    rbind(subset(dat_l_f,(cx %in% c(1,2) & lbday>=2)),
          subset(dat_l_f,cx==9) %>% mutate(lbday=2) ##,
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=7),
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=30),
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=90),
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=180),
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=270),
          ## subset(dat_l_f,cx==9) %>% mutate(lbday=360)
          ) %>% ungroup


    ## baseline flag (0 = cases, 1 = contacts, 2 = combined contacts and cases day 2)


    my_d <- augmented_dat_l_f %>%
    mutate(participant=ifelse(!cx %in% c(1,2),'contact','case')) %>%
    filter(assay %in% assays)

    my_d <- bind_rows(
        my_d %>% mutate( comp_type = participant ),
        my_d %>% filter( participant == 'contact' ) %>% mutate( comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 7, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 30, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 90, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 180, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 270, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 360, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 540, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 720, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'contact' ) %>% mutate(lbday = 900, comp_type = "'baseline'" ),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 7, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 30, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 90, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 180, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 270, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 360, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 540, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 720, comp_type = "'baseline'"),
        my_d %>% filter( participant == 'case' & lbday == 2) %>% mutate(lbday = 900, comp_type = "'baseline'"),
        ) %>%
    filter(!lbday %in% exc_lbdays)

    ## test baseline data
    for (a in seq_along(assays)){

        print(paste0("KS test for ",assays[a],":"))


        ks_test_dat <- rbind(subset(dat_l_f,(cx %in% c(1,2) & lbday>=2)),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=2),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=7),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=30),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=90),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=180),
                             subset(dat_l_f,cx==9 & lbday<=2) %>% mutate(lbday=360)
                             ) %>% ungroup %>%
                       mutate(participant=ifelse(!cx %in% c(1,2),'contact','case')) %>%
                       filter(assay %in% assays[a],lbday==2) %>%
                       select(cx,value)

        ks.test(x=ks_test_dat %>% filter(cx %in% c(1,2)) %>% ungroup %>% select(value) %>% unlist,
                y=ks_test_dat %>% filter(cx==9) %>% ungroup %>% select(value) %>% unlist) %>% print

    }

    my_d <- my_d %>% mutate(comp_type = factor(comp_type,levels=c("case","contact","'baseline'")))
    gg <- my_d %>%
    ggplot(aes(x=value)) +
    geom_density(data = my_d %>% filter(comp_type != "'baseline'"),
                 aes(fill=comp_type),col=0,bw='nrd')+
    geom_density(data = my_d %>% filter(comp_type == "'baseline'"),
                 aes(fill = comp_type),lty=2,col=1,lwd=.1,bw='nrd')+
    theme(axis.text.x = element_text(angle=90, vjust=1)) +
    labs(colour = "") +
    ylab("Density") +
    geom_rug(data=my_d  %>% filter(comp_type != "'baseline'" & participant=="case"),aes(x=value,color=comp_type),sides="t") +
    ##    geom_rug(data=my_d  %>% filter(comp_type != "'baseline'",participant=="contact"),aes(x=value,color=comp_type),sides="t") +
    facet_grid(lbday~assay,labeller = as_labeller(assay_names)) + #,scales="free_y") +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.background = element_rect(fill = "#ffffffaa", colour = NA),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          axis.text.x=element_text(size = 14),
          legend.text=element_text(size = 14)
          ) +
    scale_fill_manual(values=c(alpha("#f7f7f7",0),alpha("#fc8d62",.8), alpha("#66c2a5",.4))) +
    scale_color_manual(values=c(alpha("#fc8d62",.8), alpha("#66c2a5",.4),alpha("#f7f7f7",0))) +
    guides(color=FALSE,fill=guide_legend(override.aes=list(lty=c(2,0,0))))

    return(gg)
}



## generates and saves the full population seroprofile
## for each endmiec year at a
## specific incidence level
##' @param inc
##' @param N
##' @param endemic_yrs_per_inc
##' @return
generate_endemic_pop_profiles <- function(inc,
                                          N,
                                          endemic_yrs_per_inc){

    time_since_infection <- sapply(rep(ceiling(N*inc),endemic_yrs_per_inc), ## number of infected people
                                   function(x){
                                       last_year_inf_times <- runif(x,5,365) %>% ceiling
                                       ## assuming constant incidence historically so people infeted from
                                       ## years 1 through 2.5 in past are 1.5 times number of cases in this year
                                       inf_times_1_to_2.5_yrs <- runif(x*1.5,365,900) %>% ceiling
                                       inf_times <- c(last_year_inf_times,inf_times_1_to_2.5_yrs)
                                       c(inf_times,rep(Inf,times=N-length(inf_times)))
                                   },simplify=TRUE)
    #saveRDS(time_since_infection,paste0("generated_data/endemic_time_since_infections_full_pop_inc_",inc*1000,".rds"))

    ## generate full population sero profiles for each endemic year
    trajs <-  make_all_trajs(serotypes=c(1,2)) ## trajectories for cases

    rc <- sapply(1:ncol(time_since_infection),function(rep){
        sp <- generate_seroprof_dataset(time_since_infection[,rep],traj_tbl=trajs,round_vibs=TRUE)
        fn <- paste0("generated_data/endemic_pop_prof_inc_",1000*inc,"_",rep,".rds")
        saveRDS(sp,file = fn)
        data.frame(filename=fn,incidence = inc,stringsAsFactors = FALSE)
    },simplify=FALSE)

    return(rc)
}


## generates and saves the full population seroprofile for an epidemic
##' @title
##' @param time_since_inf
##' @param sims
##' @param inf_to_case_label
##' @return
##' @author Andrew A
generate_epidemic_pop_profiles <- function(time_since_inf,sims,inf_to_case_label){
    trajs <-  make_all_trajs(serotypes=c(1,2)) ## trajectories for cases

    rc <- sapply(1:sims,function(rep){
        generate_seroprof_dataset(time_since_inf,traj_tbl=trajs,round_vibs=TRUE) %>%
        dplyr::mutate(
                   inf_10=ifelse(t_since_inf>4 & t_since_inf<=10,1,0),
                   inf_45=ifelse(t_since_inf>4 & t_since_inf<=45,1,0),
                   inf_100=ifelse(t_since_inf>4 & t_since_inf<=100,1,0),
                   inf_200=ifelse(t_since_inf>4 & t_since_inf<=200,1,0),
                   inf_365=ifelse(t_since_inf>4 & t_since_inf<=365,1,0))
    },simplify=FALSE)

    fn <- paste0("generated_data/epidemic_pop_prof_inf_to_case_ratio_",inf_to_case_label,".rds")
    cat(sprintf("saving %s \n",fn))
    saveRDS(rc,file = fn)
    return(rc)
}

## fix mixture of log normals for contact distribution
##' @param monitor_inds - TRUE if we want to monitor the MCMC draws for each person's assignment
##' @return
fit_contact_mixture_model <- function(monitor_inds = FALSE){

    library(R2jags)
    library(runjags)

    ogawa = load_data_for_paper() %>% ungroup %>% filter(assay == "vibogaw",cx == 9,lbday2==2) %>%
    select(value) %>% unlist

    inaba = load_data_for_paper() %>% ungroup %>% filter(assay == "vibinab",cx == 9,lbday2==2) %>%
    select(value) %>% unlist

    maxvib = pmax(ogawa,inaba) %>% unname %>% log


    n_clust = 2
    clust = rep(NA,length(maxvib))
    clust[which.min(maxvib)]=1 # smallest value assigned to cluster 1
    clust[which.max(maxvib)]=2 # highest value assigned to cluster 2
    my_data = list(
        y = maxvib ,
    N = length(maxvib) ,
    n_cluster = n_clust ,
    clust = clust ,
    titers = log(c(0.01,5,10,10*2^(1:20))),
    ones_rep_Nclust = rep(1,n_clust)
    )

    inits <- function() list(p_cluster = runif(1) %>% c(.,1-.),
                             mu_of_cluster = runif(1,0,4) %>% c(.,.+2),
                             tau = runif(1,.1,1)
                             )

    set.seed(5632)
    my_monitor <-  c("p_cluster","mu_of_cluster","tau")

    ## if we want to track the chains for individual assignments
    if(monitor_inds) my_monitor <- c(my_monitor,"clust")

    my_run <- run.jags(model = "source/mixtures.jags", data = my_data,
                       inits = inits,
                       monitor=my_monitor,
                       sample = 75000,
                       n.chains = 4)

    print(summary(my_run))

    return(as.mcmc(my_run))

}


##' identify which people are in cluster 1 and which are in cluster 2
##' this re-runs the jags model
##' @author Andrew A
id_clusters <- function(){

    mcmcs <- fit_contact_mixture_model(monitor_inds=TRUE)

    clust_cols <- grep("clust\\[[0-9]+\\]",colnames(mcmcs))

    clust_assignments <- apply(mcmcs[,clust_cols] ,2,function(x) ifelse(mean(x == 1)>.5,1,2))


    rc <- bind_cols(cluster=clust_assignments,
                    id = load_data_for_paper() %>% ungroup %>% filter(assay == "vibinab",cx == 9,lbday2==2) %>% select(id))

    return(rc)
}




## helper function for keeping colors consistent through out
keep_assay_color_consistent <- function(assay_names){

    ordered(assay_names,levels = c('vibinab','vibogaw',
                                   'pglsp_s','pgctxb_s',
                                   'plpsa_s','ptctxa_s',
                                   'plpsm_s','pmctxb_s'))

    }


## makes density plot of CTB resopnses by day in
## paxvac volunteers
make_ctb_paxvac_plot <- function(){

    pv <- read_csv("data/public_challenge_data.csv")

    pv %>% filter(assay %in% c('ptctxa_s', 'pmctxb_s','pgctxb_s')) %>%
    ggplot(aes(x=log2(titer))) +
    geom_density(aes(fill=assay),alpha=.4,color=NA) +
    geom_rug() + facet_grid(day~assay,labeller = as_labeller(assay_names)) +
    theme_minimal() + guides(fill=FALSE) +
    ggtitle("North American Volunteer Post-Challenge Anti-CTB Responses")

}


##' Used to generate table with RF model performance by time period and
##' restricted to O-negatives only
##' @param igm igm subset?
##' @param target_fn filename for output
##' @return
##' @author Andrew A
make_aucs_oneg_only <- function(igm = FALSE,
                                target_fn = "table_2_oneg_only.docx"){


    if(igm){
        dat <- load_data_for_paper(long=FALSE,include_igm=TRUE) %>%
        filter(!is.na(plpsm_s) & !is.na(pmctxb_s), o_group == 0)

        full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                          'pglsp_s','pgctxb_s','pmctxb_s','plpsm_s',
                          'sex','wt','age',sep="+")


    } else {
        dat <- load_data_for_paper(long=FALSE) %>% filter(o_group == 0)

        full_set <- paste('vibinab','vibogaw','plpsa_s','ptctxa_s',
                          'pglsp_s','pgctxb_s','sex','wt','age',sep="+")

    }


    ntrees <- 1000
    k <- 20
    times <- c('inf_10','inf_45','inf_100','inf_200','inf_365')

    ## --------------------- ##
    ## First for full models ##
    ## --------------------- ##

    set.seed(262340)
    cv_roc_rfs <- lapply(1:length(times),
                         function(t){
                             my_formula <- paste0("factor(",times[t],") ~",full_set) %>% as.formula
                             print(my_formula)
                             cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)
                         })


    ## --------------------------------- ##
    ## now for reduced 2 biomarker model ##
    ## --------------------------------- ##

    set.seed(262341)
    cv_roc_rfs_reduced <- lapply(1:length(times),
                                 function(t){

                                     my_set <- paste('vibogaw','pgctxb_s',sep="+")

                                     if(t==1){
                                         my_set <- paste('vibogaw','ptctxa_s',sep="+")
                                     }

                                     my_formula <- paste0("factor(",times[t],") ~",my_set) %>% as.formula
                                     print(my_formula)

                                     cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)#,
                                     ##cutoff=c(reduced_cutoffs[t],1-reduced_cutoffs[t]))
                                 })


    ## ------------------------- ##
    ## Now for ELISA only models ##
    ## ------------------------- ##
    set.seed(264341)

    if(igm){

        ELISA_set <- paste('plpsa_s','ptctxa_s',
                           'pglsp_s','pgctxb_s',
                           'plpsm_s','pmctxb_s',sep="+")

    } else {

        ELISA_set <- paste('plpsa_s','ptctxa_s',
                           'pglsp_s','pgctxb_s',sep="+")

    }

    cv_roc_rfs_ELISA <- lapply(1:length(times),
                               function(t){
                                   my_formula <- paste0("factor(",times[t],") ~",ELISA_set) %>% as.formula
                                   print(my_formula)
                                   cv_rf_model(dat,my_formula,k=k,ntrees=ntrees)
                               })


    ## have a look at CV AUCs
    overview_full <- map_dfr(cv_roc_rfs,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_full) <- "full"


    ## have a look at performance of reduced model
    overview_reduced <- map_dfr(cv_roc_rfs_reduced,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_reduced) <- "reduced"

    overview_ELISA <- map_dfr(cv_roc_rfs_ELISA,1) %>%
    mutate(lab=paste0(round(cvAUC*100,1)," (",round(ci1*100,1),"-",round(ci2*100,1),")")) %>%
    select(lab) %>%
    t %>% unname %>% data.frame %>% set_names(paste0(c(10,45,100,200,365),"-days"))
    rownames(overview_ELISA) <- "ELISA"

    overview_tab <- (bind_rows(overview_full,overview_reduced,overview_ELISA)) %>% set_rownames(c("full model","two biomarker","ELISA only"))

    my_tab <- regulartable(bind_cols(data.frame("Model"=rownames(overview_tab)),overview_tab)) %>%
    add_header(Model = "",
               `10-days`="Infection Time Window",
               `45-days`="Infection Time Window",
               `100-days`="Infection Time Window",
               `200-days`="Infection Time Window",
               `365-days`="Infection Time Window"
               ) %>% merge_h(part = "header") %>%
    align(align="center",part="all") %>%
    autofit

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)

    print(my_tab)

}

##' Estimates log2 differences of baseline titers
##' @title
##' @return
##' @author Andrew A
estimate_baseline_titer_diffs <- function(target_fn = "tables/baseline_titer_diffs.docx"){

    dat_l_f <- load_data_for_paper(long=TRUE,full=TRUE,include_igm=TRUE) %>% mutate(u5 = ifelse(age<5,1,0) %>% factor)

    key_assays <- c('vibinab','vibogaw',
                    'pglsp_s','plpsa_s',
                    'pgctxb_s','ptctxa_s','plpsm_s')

    u5_diffs <- dat_l_f %>%
    filter(lbday==2 & assay %in% key_assays & cx %in% c(1,2)) %>% mutate(assay=clean_biomarker_labels(assay)) %>%
    group_by(assay) %>% dplyr::mutate(l2v=log2(value)) %>%
    do(w=pairwiseCI(l2v ~ u5, data=data.frame(.),method="Param.diff")) %>%
    dplyr::summarize(assay,est=w$byout[[1]]$estimate,lower=round(w$byout[[1]]$lower,2),upper=round(w$byout[[1]]$upper,2)) %>%
    mutate(`under 5` = sprintf("%.2f (%.2f, %.2f)",est,lower,upper))

    ogroup_diffs <- dat_l_f %>%
    filter(lbday==2 & assay %in% key_assays & cx %in% c(1,2)) %>% mutate(assay=clean_biomarker_labels(assay)) %>%
    group_by(assay) %>% dplyr::mutate(l2v=log2(value),o_group=factor(o_group)) %>%
    do(w=pairwiseCI(l2v ~ o_group, data=data.frame(.),method="Median.diff")) %>%
    dplyr::summarize(assay,est=w$byout[[1]]$estimate,lower=round(w$byout[[1]]$lower,2),upper=round(w$byout[[1]]$upper,2)) %>%
    mutate(`O blood group` = sprintf("%.2f (%.2f, %.2f)",est,lower,upper))

    sex_diffs <- dat_l_f %>%
    filter(lbday==2 & assay %in% key_assays & cx %in% c(1,2)) %>% mutate(assay=clean_biomarker_labels(assay)) %>%
    group_by(assay) %>% dplyr::mutate(l2v=log2(value),sex=factor(sex)) %>%
    do(w=pairwiseCI(l2v ~ sex, data=data.frame(.),method="Median.diff")) %>%
    dplyr::summarize(assay,est=w$byout[[1]]$estimate,lower=round(w$byout[[1]]$lower,2),upper=round(w$byout[[1]]$upper,2)) %>%
    mutate(sex = sprintf("%.2f (%.2f, %.2f)",est,lower,upper))

    my_tab <- full_join(
        u5_diffs %>% select(assay,`under 5`),
        ogroup_diffs %>% select(assay,`O blood group`)) %>%
    full_join(sex_diffs %>% select(assay,sex)) %>% regulartable %>% autofit()

    doc <- read_docx()
    doc <- body_add_flextable(doc, value = my_tab)
    print(doc, target = target_fn)


}
