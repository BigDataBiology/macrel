<template>
    <div class="amps-container">
        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="2"><div class="grid-content"></div></el-col>
            <el-col :span="14">
                <div class="grid-content amps-head">
                    <h3>result of AMP prediction</h3>
                    <el-divider></el-divider>

                    <div class="resultArea">
                        <div class="panel-head">
                            <h3>classification</h3>
                            <el-button type="primary" @click="downloadResult">download<i class="el-icon-download el-icon--right"></i></el-button>
<!--                            <el-button type="success" icon="el-icon-download" circle @click="downloadResult"></el-button>-->
                        </div>
                        <div class="panel-body">
                            <el-table
                                    v-loading="loading"
                                    element-loading-text="loading"
                                    element-loading-spinner="el-icon-loading"
                                    element-loading-background="rgba(0, 0, 0, 0.5)"
                                    :data="ampList"
                                    style="width: 100%"
                                    max-height="300"
                                    stripe
                                    empty-text = "empty"
                                    :default-sort = "{prop: 'access', order: 'descending'}"
                            >
                                <el-table-column type="index" :index="indexMethod"></el-table-column>
                                <el-table-column
                                        sortable
                                        prop="access"
                                        label="access">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        prop="amp_family"
                                        label="amp_family">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        prop="amp_probability"
                                        label="amp_probability">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        prop="hemolytic"
                                        label="hemolytic">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        prop="hemolytic_probability"
                                        label="hemolytic_probability">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        prop="sequence"
                                        label="sequence">
                                </el-table-column>
                            </el-table>
                        </div>
                        <div class="panel-tail">
                            <h3>{{ampList.length}}&nbsp;AMPs</h3>
                        </div>
                    </div>
                </div>
            </el-col>
            <el-col :span="2"><div class="grid-content"></div></el-col>

        </el-row>
    </div>
</template>

<script>
    export default {
        name: "AMPs",

        data(){
            return {
                filePath:'',
                ampList:[],
                loading: true,
            }
        },

        // 页面创建时，要进行的操作
        created() {
            this.loadResultObject();
            this.loading = false;
        },

        methods:{
            async loadResultObject(){
                var resultObjectStr = await window.sessionStorage.getItem('resultObjectStr');
                let resultObject = JSON.parse(resultObjectStr);
                if (resultObject.code !== 1){
                    return this.$message.error("load data error.");
                }

                this.ampList = resultObject.data.objects;
                this.filePath = resultObject.data.filePath;
            },

            indexMethod(index) {
                return index + 1;
            },

            async downloadResult(){
                if (this.ampList.length<=0){
                    this.$message.error("no data!");
                    return;
                }

                await this.$http({
                    method:'get',
                    url:'/file/download',
                    params:{
                        filePath:this.filePath,
                    },
                    headers:{
                        'Content-Type':'application/x-download;charset=utf-8',
                    },
                    responseType:'blob',
                })
                .then( response => {
                    this.download(response);
                })
                .catch( (e)=> {
                    this.$message.error('Access failed');
                    window.console.log(e);
                } );
            },

            async download(response){
                if (!response) {
                    return this.$message.error('file not exits.');
                }
                let url = window.URL.createObjectURL(new Blob([response.data]));
                let link = document.createElement('a');
                link.style.display = 'none';
                link.href = url;
                link.setAttribute('download',response.headers['filename']);

                document.body.appendChild(link);
                link.click();

                document.body.removeChild(link); //下载完成移除元素
                window.URL.revokeObjectURL(url); //释放掉blob对象
            },

        },

    }
</script>

<style scoped lang="less">

    .amps-container{
        margin-top: 0.1em;
        height: 100%;
    }

    .resultArea{
        border: 1px solid;
    }

    .panel-head{
        justify-content: space-between;
        display: flex;
        align-items: center;
        padding: 0 6em;
        height: 4em;
        background-color: antiquewhite;
    }

    .panel-body{
        padding: 1em 2em;
    }

    .panel-tail{
        display: flex;
        height: 2em;
        align-items: center;
        padding: 0 10em;
        margin: 2em 0;
        justify-content: flex-end;
    }

    .el-row {
        height: 100%;
        padding: 1em 0;
        /*background-color: greenyellow;*/
    }

    .el-table{
        margin-top: 1em;
        font-size: 1.1em;
    }

</style>