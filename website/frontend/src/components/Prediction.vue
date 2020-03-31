<template>

    <div class="prediction_container">

        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="1"><div class="grid-content"></div></el-col>
            <el-col :span="18" class="el-col-form">
                <div class="grid-content formArea">
                    <div class="form-head">Antimicrobial activity prediction</div>
                    <div class="form-body">
                        <el-form ref="pFormRef" :model="pFormModel" :rules="pFormRules" size="medium" class="prediction_form">
                            <div style="float: left; width: 60%">
                            <h4 style="margin-top: 0.4em;margin-bottom: 0.2em">You can paste peptides or contigs sequence in <span style="color: red">FASTA</span> format into the field below:</h4>

                            <el-form-item prop="textData" style="margin-bottom: 0.1em">
                                <el-input type="textarea" :rows="8" v-model="pFormModel.textData" placeholder="input your sequence here"></el-input>
                            </el-form-item>
                            <p style="text-align: right; margin-top: 0px; padding-top: 0px; font-size: small"><a href="#" @click="doExample('peptides')">Peptides example</a> || <a href="#" @click="doExample('contigs')">Contigs example</a></p>

                            </div>
                            <div style="float: right; width: 35%">
                            <h4 style="margin-top: 0.4em;margin-bottom: 0.2em">Or submit a file in FASTA format:</h4>

                            <el-form-item>
                                <!--
                                A colon is added to indicate that a variable or expression follows; without a colon is a corresponding string literal!
                                :auto-upload=false  // No automatic upload
                                :limit=1
                                drag
                                action="" // Do not set the upload address for the time being, because we want to intercept the upload in the form
                                -->
                                <el-upload
                                        ref="upload"
                                        action=""
                                        show-file-list
                                        class="file-upload"
                                        drag
                                        :auto-upload="false"
                                        accept=".fasta,.fa,.faa,.fna"
                                        :on-exceed="handleExceed"
                                        :on-remove="handleRemove"
                                        :multiple="false"
                                        :limit="1"
                                        >
                                    <i class="el-icon-upload"></i>
                                    <div class="el-upload__text">Drag or click here to upload</div>
                                    <div class="el-upload__tip">(FASTA format, limited to 500KB)</div>
                                </el-upload>
                            </el-form-item>
                            </div>

                            <div style="clear: both" />

                            <el-divider><i class="el-icon-more-outline"></i></el-divider>
                            <el-form-item label="Data Type: " required prop="dataType" style="margin-bottom: 0.3em">
                                <el-radio-group v-model="pFormModel.dataType">
                                    <el-radio label="peptides">Peptides</el-radio>
                                    <el-radio label="contigs">Contigs (nucleotide)</el-radio>
                                </el-radio-group>
                            </el-form-item>

                            <br/>
                            <el-form-item class="btns">
                                <el-button type="primary" @click="onSubmit('pFormRef')">Submit</el-button>
                                <el-button type="danger" @click="resetForm('pFormRef')">Clear Form</el-button>
                            </el-form-item>
                        </el-form>
                    </div>

                    <p>For large inputs (>500KB or >1,000 sequences), please <a
                    href="https://github.com/BigDataBiology/macrel/">download
                    the tool</a> and run it locally.</p>

                </div>
                <div class="introducion">
                    <p>Antimicrobial peptides (AMPs) are small proteins (10–100
                    amino acids) that can either cause lysis or generally
                    interfere with cell growth.The rapid growth of
                    antibiotic-resistant microorganisms in the last years
                    became a world health problem. Thus, AMPs could represent a
                    source of new antibiotic treatments. Recently, more and
                    more genomic and metagenomic sequences have become publicly
                    available. Macrel is pipeline that overcomes most of the
                    problems associated to the prospection of AMPs in high
                    throughput data sets analyzing diverse types of inputs
                    (reads, contigs, protein sequences). The protein sequences
                    predicted as AMPs by Macrel are, usually, specific and
                    accurate, what can reduce the cost of experiments.</p>

                    <p>Macrel uses machine learning to select peptides with
                    high probability of being an AMP. Furthermore, Macrel is
                    also capable to perform the classification of AMPs into
                    hemolytic and non-hemolytic peptides. This allows
                    researchers to select the most interesting peptides for
                    further testing <i>in vitro</i>.</p>

                    <p>Peptides submitted to the Macrel prediction should
                    consist of 20 canonical amino acids and their length should
                    range from 10 to 100 amino acids. Please avoid contigs
                    containing non-canonical bases, such as N, R or Y.</p>

                    <p>This prediction system will output a table containing
                    the predicted AMP sequence, its probability to be an AMP
                    and the family it belongs, as well as the predictions
                    relative to the hemolytic activity.</p>

                </div>
                <div class="glossary">
                    <h3>Glossary:</h3>
                    <ul>
                        <li>ALP - Anionic linear peptides</li>
                        <li>ADP - Anionic dissulphide-bond forming peptides</li>
                        <li>CLP - Cationic linear peptides</li>
                        <li>CDP - Cationic dissulphide-bond forming peptides</li>
                        <li>HEMO – Hemolytic peptide</li>
                        <li>NonHEMO – Non-Hemolytic peptide</li>
                    </ul>
                </div>
            </el-col>
            <el-col :span="1"><div class="grid-content"></div></el-col>
        </el-row>
    </div>

</template>

<script>
    export default {
        name: "Prediction",

        data() {
            return {

                pFormModel: {
                    textData: '',
                    dataType: '',
                },

                pFormRules: {

                    dataType: [
                        { required: true, message: 'Please select data type', trigger: 'change' },
                    ],

                },

            };
        },

        methods: {

            handleRemove(){
                console.log("handleRemove");
                //refresh the browser to clear file cache.
                this.$router.go(0);
            },

            onSubmit(formRef) {

                this.$refs[formRef].validate(async valid => {

                    if (!valid) return false;

                    let pForm = this.$refs[formRef].$el;
                    let formData = new FormData(pForm);
                    //console.log(formData.get("file"));

                    //check whether textData and file both are null.
                    let fileSize=formData.get("file").size;
                    if (!this.pFormModel.textData && fileSize===0) {
                        this.$message.error({message: 'Sequence is empty and no file uploaded!', duration:5000});
                        return false;
                    }

                    //limit size of uploaded file
                    if(fileSize>0){
                        const isLimit = fileSize / 1024  < 500;
                        if (!isLimit) {
                            this.$message.error({message: 'File size exceeds 500KB. Please download the tool and run locally for large inputs', duration:5000});
                            return false;
                        }
                    }

                    //limit textData rows
                    if(this.pFormModel.textData){
                        var array=this.pFormModel.textData.split("\n");

                        var sequenceLength=array.length;
                        if (sequenceLength > 1000){
                            this.$message.error({message: 'Sequence exceeds 1000 rows. Please download the tool and run locally for large inputs', duration:3000})
                            return false;
                        }
                    }

                    formData.append('textData', this.pFormModel.textData);
                    formData.append('dataType', this.pFormModel.dataType);

                    let config = {
                        headers: {
                            'Content-Type': 'multipart/form-data'
                        }
                    };

                    const loading = this.$loading({
                        lock: true,
                        text: 'Analyzing, please wait',
                        spinner: 'el-icon-loading',
                        background: 'rgba(0, 0, 0, 0.5)',
                    });

                    let resultObject; // This will hold the result
                    await this.$http.post('http://localhost:5000/predict', formData, config)
                        .then(
                            response => {
                                if (response.status === 200){
                                    resultObject = response.data;
                                    window.console.log("submit success");
                                }
                            }
                        ).catch( ()=> {
                            loading.close();
                            this.$message.error({message: 'No response from API server. Please try again later',duration:5000});
                        })
                    ;

                    loading.close();
                    if (resultObject.code !==1){
                        this.$message.error(resultObject.message);
                        return ;
                    }

                    // Turn the resulting object into a string and save it in a browser session
                    let resultObjectStr = JSON.stringify(resultObject);
                    window.sessionStorage.setItem('resultObjectStr',resultObjectStr);

                    await this.$router.push({name: "amps"});
                });
            },

            handleExceed() {
                this.$message.warning({message: 'Only one file can be uploaded at a time', duration:5000});
            },

            resetForm(formRef) {

                let pForm = this.$refs[formRef].$el;
                let formData = new FormData(pForm);

                if(formData.get("file").size>0){
                    //refresh the browser to clear file cache.
                    this.$router.go(0);
                }else{
                    this.$refs['upload'].clearFiles();
                    this.$refs[formRef].resetFields();
                }
            },

            doExample(etype) {
                window.console.log("doExample");
                this.pFormModel.dataType = etype;
                this.onExample('pFormRef');
            },

            onExample(refName){
                this.$refs[refName].validate(valid => {
                    if(!valid) return;
                    var dataType = this.pFormModel.dataType;
                    // console.log(dataType);
                    if (dataType === 'peptides'){
                        this.pFormModel.textData = ">AP00002|AMP\n" +
                            "YVPLPNVPQPGRRPFPTFPGQGPFNPKIKWPQGY\n" +
                            ">AP00007|AMP\n" +
                            "GNNRPVYIPQPRPPHPRL\n" +
                            ">P20491|NAMP\n" +
                            "MISAVILFLLLLVEQAAALGEPQLCYILDAVLFLYGIVLTLLYCRLKIQVRKAAIASREKADAVYTGLNTRSQETYETLKHEKPPQ\n" +
                            ">P35109|NAMP\n" +
                            "MVDDPNKVWPTGLTIAESEELHKHVIDGSRIFVAIAIVAHFLAYVYSPWLH\n" +
                            ">P19962|NAMP\n" +
                            "DSDSAQNLIG";
                    }
                    if (dataType === 'contigs'){
                        this.pFormModel.textData = '>scaffold2530_2_MH0058\n' +
                            'CTTCTGATCTTTACGCAGCATTGTGTGTTTCCACCTTTCAAAAAATTCTCCGTGAACTGC\n' +
                            'GCCCTGGGAGTGGTGAAATCCTCCGCGGAACGAAGTCCCGGAATTGCGCACAAATTCACG\n' +
                            'TGCTGAACAATTTTACCATAGGAATGTGCGGTTGTAAAGAGAAAAATGCAAAAAATTCCT\n' +
                            'TATTTTTATAAAAGGAGCGGGGAAAAGAGGCGGAAAATATTTTTTTGAAAGGGGATTGAC\n' +
                            'AGAGAGAAACGGCCGTGTTATCCTAACTGTACTAACACACATAGTACAGTTGGTACAGTT\n' +
                            'CGGAGGAACGTTATGAAGGTCATCAAGAAGGTAGTAGCCGCCCTGATGGTGCTGGGAGCA\n' +
                            'CTGGCGGCGCTGACGGTAGGCGTGGTTTTGAAGCCGGGCCGGAAAGGAGACGAAACATGA\n' +
                            'TGCTGTTTGGTTTTGCGGGGATCGCCGCCATCGTGGGTCTGATTTTTGCCGCTGTTGTTC\n' +
                            'TGGTGTCCGTGGCCTTGCAGCCCTGAGAACGGGGCAGATGCAATGAGTACGCTGTTTTTG\n' +
                            'CTTGGTATCGCAGGCGCGGTACTGCTGGTTATTTTGCTGACAGCGGTGATCCTGCACCGC\n' +
                            'TGATCGAACATTCCTCAGAAAGGAGAGGCACACGTTCTGACATTGAATTACCGGGATTCC\n' +
                            'CGTCCCATTTATGAACAGATCAAGGACGGCCTGCGGCGGATGATCGTCACCGGGGCC\n' +
                            '>scaffold75334_1_MH0058\n' +
                            'ACAGCTTCCTCACCATCAACAGCCACTGCTACGATACCGGCAGGAACAGAGATTGTAGCG\n' +
                            'TTATCGGAAGTAAGAACGGTCTCAGCGATTACCTTACCCAAGATATTCGTGATAACTACA\n' +
                            'GACTTACCAGCTGCGCCTTGAACGGTTACGGTACCGTTACCGGCTACTACAGAGATACCT\n' +
                            'TCTACAGCATCGATTGTTTCGTTATCTGTTGCGATATCATCACCTTGCTCCACATTAAAG\n' +
                            'ATCAAAGCATCGTCACCACCTGTATTCATCTCATCGAAATTAGAAGAACGATCATCAGAC\n' +
                            'AGAACTAAACAACCATTCTGCATTCTCAACCAAGCCGCATATGTAGGAGCGATGTCACCC\n' +
                            'ATAGCAGAATTATTGTAACCAGGAACATGTTTACTCTTCATAGACTCGAACAAGAATGCT\n' +
                            'CTATCAGCTTCTACCTCATTAGCGGCAACTTCCTTGTTTACGAAACGCATAGACCAAGTC\n' +
                            'ACATACTTATGGTTATCACCAGATAAGATATATTTGTGTGAAACAAATTCATCTTCAGAT\n' +
                            'ACTCCAGCCTTCTTAGCCGCAGTCCAAGCATCCTCCTCAGCCTTGTTCAAATCAGCGAAG\n' +
                            'TTGATCTTCTCGTTCGGTAAGTTCTTGAACTCATCTCTCAAGATATACAAAGTATCAGCT\n' +
                            'ACACGGATAGCTTCCACGAAACCTGCACGATCATATTTCTTCCACATGTAGTCAGCATCT\n' +
                            'TGCTCCCCATTTGATAATGTTGACCAACCACATTTATTAGCGAAATCATGGAAATTGATC\n' +
                            'AAATATTTACCTCTTTCAAAGCCCGGAACAGCCGGTTTTGCGTGAACGCAATGCCATTTA\n' +
                            'TCTGTCGGTTTACCTTCTGCGGTAAAGTGATGGTCATCTTCCGTACAAGGTACGCCTTCA\n' +
                            'ACTCCTTCGAAATCGTTACGATCAATAGAGATCAAATATTGAGGCTTGATATTACCAGAA\n' +
                            'CCACGGTTTACCCATGCGGTATCAACGATGAATGACAAACCATCCTCAGTCTTATCCGGA\n' +
                            'GTATAAATACCTAAGAAGTCAATACCTTCTTTCATGAAGTTCTTATTATTCTCAACCTGC\n' +
                            'AAGTATTCTTTACGGTATTTCTCGATAAAACGTAAAGTATCGGCCTTGTCACCTTCATTA\n' +
                            'CCTTCAAGTTCAAGAGAACTGAAACGACGGTAAAGCGGAGTGTTGTCCGGTTCGATAGCG\n' +
                            'AAAGCAGAAGTACGAGTCTCGCCTAATACTTGATTCTTCAATGTAGCGGCACCGTCATAA\n' +
                            'TCGGATACTCCAGCTTTCTGATAACCTAATACATATTTTGTAGCATTAGAATAATTATCA\n' +
                            'TATGCAGCATTAACAATCGCATAATAGTGCTTTCCACTGATATAGTTGTTTTCTTTGAAG\n' +
                            'AATGCATAATAAGTCGTATAAGATTGCTGATGACCAGTAGCACTAGAACCTACATTAAAC\n' +
                            'TTATCTTCCTTATTAACTTTCCATTCATTAAGTTCACCTTTACCATCTTTATAAGAAACC\n' +
                            'TCATAAGCAGTTCTCACTAGAGTCTTCAAGCCTTTAATACGATTCTTAGCCACATCATCG\n' +
                            'CTAACCTTATAACCGTAAGGGATATCTACAGTCGCAAATATCTTTGAGTTATTATTCCAA\n' +
                            'CGCTCAGTCAATGTATCAATTACGAAAGCATCTTTTCCATCCAACACAGTCAACGTAGAG\n' +
                            'TCTGTAGATTTAGCGATGTATTTATCCGTAGCATAAGGATGCCAATAGTTAAACGTATAT\n' +
                            'TTGTTAACCAACAAAGAATCTTTCTCCAACTTTCTGTAACCTAAATACGGATCAGAAATA\n' +
                            'GAAGCTGCCGGAACTTTCTCGAAGAACAAAGAATCTATACCGACAGAGCCATACTCCGCT\n' +
                            'GTAAGTTTATCTCCATATACAAACATATAAGAACCACCCTCTGCCTTTCTTAACTGAACA\n' +
                            'GTAGAATAAACAGCAGTCTTATATGACTCCATCTGGCCATCACCATCAAAATCCGCATTT\n' +
                            'ATAACTTTATCTGCGAATTCACGGTTAGAGATAGTAACCGGAGACACAGCCTTCACTTTA\n' +
                            'TCATTACTTGTATTTTTCTTCAATACAACCCATTGATAAGCTGGCATATGAGCAACACTC\n' +
                            'TGCTCATCCTCATTCACGGTTGTCCAACGGATAGTACCATTCTCATAGATAGGAGAAGCT\n' +
                            'AAATACTGTCCTTGTTTATTCTTGATGAAATAAACACCATCATCAACGGAAGTCTTGGAA\n' +
                            'GAACCCTTTTCCTCACATCCGCTAAGACCCAAAGAAATATGAGTGTTCTGGTCCTCCGGA\n' +
                            'TCAACTGTCAAAATACAAGTTTCATCTTTAATCAAGTCTTGCAAAGAGACACGCCAGAAC\n' +
                            'TCATGCTCGTTAACTGGAGTAGCCGTAACGTGCGTATTCCAATACTTGCCATCATTATTC\n' +
                            'CAAGTAACTTTTTTTACATAGATGTATAAGCTATCGCCACTCGGAGAATAGATGAATTTG\n' +
                            'AACTGATGTTGAGCACGTAATTCAGAGCTAATAGAACTGTAAGCATCAGCTTTATCCTTT\n' +
                            'GCTTTCTCTGTCCAACCGTAAGCTAAGAACTTTGTACCTGTCTCATTCGTATAAGCCGTA\n' +
                            'TCAACTTTCAAATAAGCGTCTTTATCCTTTTGCTTGATGAAAAGCCATTTGTTATTATCT\n' +
                            'TTATCCTCTGCGATAAACTTCTGCTCATTGAACGGATTCTTTAAGCTAGTACCTTTCACG\n' +
                            'TCCGGAGTGAAAGTCAATTGTACTCCAGCCTTGTTATCTTGAAGAAGACCCAACTTCGTA\n' +
                            'TTTACTTGATCTTCGCTCAATGCGATTTGTGCTGCGCCAACCAAAGCAAAGTTAGTCACA\n' +
                            'TCAGCGGCGTTTGAAGTGGCATCTATTTGATTAGCACCCCACTTGGCAACCTTAACAGCT\n' +
                            'CCGGTAGTAGGATCCGCTTTCAAACCAACAATTGAGTCCGTTGAGAAATAAGAATACAAA\n' +
                            'GGTCTCTTTTCTTCCAAATCCTTGTATGAGCGAGAGAATGCCCAACCTGCGATCTCACCA\n' +
                            'CCTAAAACAGGTTTCCAAGCACCATTCGTACCATTCTCAGTCTTCTCATGGCCAGCCATT\n' +
                            'GTCAAATCCAACATAGTTCCAGACAATTTATTCTGGAAGTCATAAATAGGAGCTTGACCT\n' +
                            'TGATTATAATTAGAGACAGAAACACACCACAATGTAGCCTCTAAGTTATTTTTAGCTTCG\n' +
                            'CTTGCACTATAAACACGTAGTTCATAATCACCTGTTCCACGATTTAGCTCCATCGCAAGA\n' +
                            'TAAGCAGGAGTAGTGCCATCCATGACCTTCAACTGATAAAGACCAGAATTAGCTCCCTCC\n' +
                            'TTCAAACCGCCTAAACGCCATTCTTGACCCAATACAGTTTCTGGGTCTACTCGACTTGGT\n' +
                            'GTAGTCTGTGCATTAACAGACATAACACTTAACAGTGCCATACCTGCCAAAAGAGTAGAA\n' +
                            'AACT';
                    }
                });
            }

        }
    }
</script>

<style scoped lang="less">
    .prediction_container{
        height: 100%;
    }

    .el-row {
        height: 100%;
        padding: 1rem 0;
        /*background-color: #f9fafc;*/
    }

    .formArea{
        border: 1px solid;
    }

    .form-head{
        border-bottom: 1px solid;
        background-color: antiquewhite;
        padding: 1em 2em;
    }

    .form-body{
        padding: 1em 0;
    }

    .el-upload-dragger{
        width: 200px;
    }

    .prediction_form{
        padding: 0 2em;
        box-sizing: border-box;
    }

    .btns{
        display: flex;
    }

    .introducion{
        margin-top: 2em;
    }

    p{
        line-height: 1.5em;
        text-indent: 3em;
    }

    .glossary{
        margin: 2em 0;
        padding: 1em;
        background-color: #f8f9fa;
    }

    li{
        line-height: 2em;
    }


</style>
