<template>

    <div class="choice_container">

        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="1"><div class="grid-content"></div></el-col>
            <el-col :span="18">
                <div class="cho">

                    <router-link :to="{name:'prediction',params:{dataType:'contigs'}}">
                        <el-button round>
                            <img title="DNA" src="../assets/image/icons8-dna-96.png"/>
                        </el-button>
                    </router-link>

                    <router-link :to="{name:'prediction',params:{dataType:'peptides'}}">
                        <el-button round>
                            <img title="Peptides" src="../assets/image/icons8-protein-96.png"/>
                        </el-button>
                    </router-link>

                </div>

            </el-col>
            <el-col :span="1"><div class="grid-content"></div></el-col>
        </el-row>
    </div>

</template>

<script>
    export default {
        name: "Choice",

        data() {
            return {

            };
        },

        methods: {

        }
    }
</script>

<style scoped lang="less">
    .choice_container{
        margin-top: 2em;
    }
    .cho{
        display: flex;
        justify-content: space-evenly;
    }
</style>
