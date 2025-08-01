#!/usr/bin/env nextflow

/*
    CDC Tick Surveillance Amplicon Sequencing Analysis Pipeline

    December 13, 2021

    Mark Stenglein
*/

nextflow.enable.dsl = 2                                                         

// use nf-core style initialization and checking of workflow parameters
WorkflowMain.initialise(workflow, params, log)                                  

/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE                                                 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
include { MPAS_WORKFLOW } from './workflows/mpas'                                      
                                                                                
/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS                                                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
//                                                                              
// WORKFLOW: Execute a single named workflow for the pipeline                   
// See: https://github.com/nf-core/rnaseq/issues/619                            
//                                                                              
workflow {                                                                      
    MPAS_WORKFLOW ()                                                       
}     

