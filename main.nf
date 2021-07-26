#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/conva
========================================================================================
 nf-core/conva Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/kaurravneet4123
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

log.info Utils.logo(workflow, params.monochrome_logs)

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/conva --input samplesheet.csv -profile conda"
    log.info NfcoreSchema.paramsHelp(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --        FASTA PARAMETER VALUE             -- */
////////////////////////////////////////////////////

params.fasta = Workflow.getGenomeAttribute(params, 'fasta')



////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params, json_schema)
log.info NfcoreSchema.paramsSummaryLog(workflow, params, json_schema)
log.info Workflow.citation(workflow)
log.info Utils.dashedLine(params.monochrome_logs)

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////

Workflow.validateMainParams(workflow, params, json_schema, log)

////////////////////////////////////////////////////
/* --            RUN WORKFLOW(S)               -- */
////////////////////////////////////////////////////
include { CONVA } from './workflows/conva' addParams( summary_params: summary_params )

workflow {
    CONVA ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
