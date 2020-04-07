<template>
    <div class="amps-container">
        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="1"><div class="grid-content"></div></el-col>
            <el-col :span="18">
                <div class="grid-content amps-head">
                    <h3>Result of AMP prediction</h3>
                    <p><strong>Note that only predicted AMPs are displayed.</strong></p>
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
                                    :data="currentPageData"
                                    style="width: 100%"
                                    max-height="330"
                                    stripe
                                    fit
                                    empty-text = "empty"
                                    :default-sort = "{prop: 'Access', order: 'descending'}"
                            >
                                <el-table-column type="index" :index="indexMethod" width="60" header-align="center" align="center">

                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        align="center"
                                        height="70"
                                        prop="Access"
                                        width="150"
                                        label="Peptide ID">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        align="center"
                                        width="80"
                                        prop="AMP_family"
                                        label="Family">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        width="150"
                                        prop="AMP_probability"
                                        label="Probability">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        align="center"
                                        width="150"
                                        prop="Hemolytic"
                                        label="Is hemolytic">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        width="150"
                                        prop="Hemolytic_probability"
                                        label="Is hemolytic (probability)">
                                </el-table-column>
                                <el-table-column
                                        sortable
                                        header-align="center"
                                        align="center"
                                        width="200"
                                        prop="Sequence"
                                        label="sequence">
                                </el-table-column>
                            </el-table>
                            <div class="paginationClass">
                                <el-pagination
                                        background
                                        pager-count="4"
                                        hide-on-single-page
                                        @size-change="handleSizeChange"
                                        @current-change="handleCurrentChange"
                                        layout="total, sizes, prev, pager, next, jumper"
                                        :total=totalItems
                                        :page-sizes=pageSizes
                                        :page-size=pageSize
                                        :current-page=currentPage
                                >
                                </el-pagination>
                            </div>
                        </div>

                        <div class="panel-tail"></div>
                    </div>
                    <p>Results with macrel version <i>{{macrel_version}}</i></p>
                </div>
            </el-col>
            <el-col :span="1"><div class="grid-content"></div></el-col>

        </el-row>
    </div>
</template>

<script>
    export default {
        name: "AMPs",

        data(){
            return {
                totalItems:0,
                pageSizes:[10, 20, 50, 100],
                pageSize:5,
                currentPage:1,
                currentPageData:[],

                ampList:[],
                macrel_version: '?',
                loading: true,
            }
        },

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
                this.macrel_version = resultObject.macrel_version;
                this.rawdata = resultObject.rawdata;

                //first page
                this.totalItems=this.ampList.length;
                if(this.totalItems>this.pageSize){
                    for (let index = 0; index < this.pageSize; index++) {
                        this.currentPageData.push(this.ampList[index]);
                    }
                }else{
                    this.currentPageData = this.ampList;
                }
            },

            indexMethod(index) {
                return index + 1;
            },

            //change size of every page
            handleSizeChange(pageSize) {
                this.pageSize = pageSize;
                this.handleCurrentChange(this.currentPage);
            },

            //change page number
            handleCurrentChange(currentPage) {
                this.currentPage = currentPage;
                this.currentChangePage(this.ampList);

            },

            //update page
            currentChangePage(list) {
                let from = (this.currentPage - 1) * this.pageSize;
                let to = this.currentPage * this.pageSize;
                this.currentPageData = [];
                for (; from < to; from++) {
                    if (list[from]) {
                        this.currentPageData.push(list[from]);
                    }
                }
            },

            async downloadResult(){
                if (this.ampList.length<=0){
                    this.$message.error("no data!");
                    return;
                }

                // The code below simulates a download
                let url = window.URL.createObjectURL(new Blob([this.rawdata]));
                let link = document.createElement('a');
                link.style.display = 'none';
                link.href = url;
                link.setAttribute('download', 'macrel.out.tsv');

                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                window.URL.revokeObjectURL(url);
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



    .el-row {
        height: 100%;
        padding: 1em 0;
        /*background-color: greenyellow;*/
    }

    .el-table{
        margin-top: 1em;
        font-size: 1.1em;
    }

    .paginationClass{
        display: flex;
        height: 2em;
        align-items: center;
        margin: 2em 0;
        justify-content: flex-start;
    }

</style>
