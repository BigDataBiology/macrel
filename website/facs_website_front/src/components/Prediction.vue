<template>

    <div class="prediction_container">

        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="2"><div class="grid-content"></div></el-col>
            <el-col :span="14">
                <div class="grid-content formArea">
                    <div class="form-head">This is Antimicrobial peptides prediction.</div>
                    <div class="form-body">
                        <el-form ref="pFormRef" :model="pFormModel" :rules="pFormRules" size="medium" class="prediction_form">
                            <h4 style="margin-top: 0.4em;margin-bottom: 0.2em">You can paste peptides or contigs sequence(less than 500 rows) with <span style="color: red">FASTA/FA</span> format into the field below:</h4>

                            <el-form-item prop="textData" style="margin-bottom: 0.1em">
                                <el-input type="textarea" :rows="8" v-model="pFormModel.textData" placeholder="input your sequence here"></el-input>
                            </el-form-item>
                            <br/>

                            <h4 style="margin-top: 0.5em;margin-bottom: 0.2em">Or submit a file(less than 500KB) in FASTA format:</h4>

                            <el-form-item>
                                <!--
                                加冒号的，说明后面的是一个变量或者表达式；没加冒号的后面就是对应的字符串字面量！
                                :auto-upload=false  // 取消自动上传
                                :limit=1 // 限制只能上传一张，
                                drag // 设置这个让可以把图片拖进来上传
                                action="" // 暂时不设置上传地址，因为我们是要拦截在form中上传
                                -->
                                <el-upload
                                        ref="upload"
                                        class="file-upload"
                                        :auto-upload="false"
                                        show-file-list
                                        drag
                                        action=""
                                        accept=".fasta,.fa,.faa,.fna"
                                        :on-success="onsuccess"
                                        :on-exceed="handleExceed"
                                        :multiple="false"
                                        :limit="1"
                                        >
                                    <i class="el-icon-upload"></i>
                                    <div class="el-upload__text">Drag files here，or <em>Click upload</em></div>
                                    <div class="el-upload__tip" slot="tip">only FASTA/FA file ，and no more than 500KB</div>
                                </el-upload>
                            </el-form-item>

                            <el-divider><i class="el-icon-more-outline"></i></el-divider>
                            <el-form-item label="Data Type" required prop="dataType" style="margin-bottom: 0.3em">
                                <br/>
                                <el-radio-group v-model="pFormModel.dataType">
                                    <el-radio label="peptides"></el-radio>
                                    <el-radio label="contigs"></el-radio>
                                </el-radio-group>
                            </el-form-item>

                            <br/>
                            <el-form-item class="btns">
                                <el-button type="primary" @click="onSubmit('pFormRef')">Submit</el-button>
                                <el-button type="danger" @click="resetForm('pFormRef')">Clear Form</el-button>
                                <el-button type="info" @click="onExample('pFormRef')">Example</el-button>
                            </el-form-item>
                        </el-form>
                    </div>
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
            <el-col :span="2"><div class="grid-content"></div></el-col>
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

                fileList:[],

                pFormRules: {

                    dataType: [
                        { required: true, message: 'please select data type', trigger: 'change' },
                    ],
                },


            };
        },

        methods: {

            //阻止upload的自己上传，进行再操作
            // eslint-disable-next-line no-unused-vars
            onsuccess(file,fileList) {
                //重新写一个表单上传的方法
            },

            onSubmit(formRef) {

                this.$refs[formRef].validate(async valid => {
                    if (!valid) return;

                    let pForm = this.$refs[formRef].$el;
                    let formData = new FormData(pForm);
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
                            this.$message.error({message:'The server is under maintenance,please try again later',duration:5000});
                        })
                    ;

                    loading.close();
                    if (resultObject.code !==1){
                        this.$message.error(resultObject.msg);
                        return ;
                    }

                    // Turn the resulting object into a string and save it in a browser session
                    let resultObjectStr = JSON.stringify(resultObject);
                    window.sessionStorage.setItem('resultObjectStr',resultObjectStr);

                    this.$router.push({path:"/prediction/amps"});
                });
            },

            handleExceed() {
                this.$message.warning({message:'Only one file can be uploaded at a time',duration:5000});
            },


            resetForm(formRef) {
                this.$refs['upload'].clearFiles();
                this.$refs[formRef].resetFields();
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
