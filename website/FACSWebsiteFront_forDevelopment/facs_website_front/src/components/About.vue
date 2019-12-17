<template>
    <div class="about-container">
        <el-row type="flex" class="row-bg" justify="space-around">
            <el-col :span="2"><div class="grid-content"></div></el-col>
            <el-col :span="14">
                <div class="content">
                    <div class="content-header">
                        <h3 style="margin: 0">About FACS</h3>
                    </div>
                    <div class="content-body">
                        <p>The Fast Antimicrobial peptides Classification System (FACS) is developed at the Institute of Science and Technology for Brain-inspired Intelligence (ISTBI) at Fudan University, Shanghai (China).</p>
                        <p>The classifiers implemented in FACS use random forests and a short set of 22 descriptors to predict AMPs with an accuracy above 96% and 100% of precision. FACS was designed to find novel peptides in metagenomic samples, an important resource of diversity. Please visit our github webpage if you are interested in running it locally.</p>
                        <el-divider></el-divider>
                        <h3>Funding:</h3>
                        <p>This project was supported by Shanghai Municipal Science and Technology Major Project (2018SHZDZX01), ZHANGJIANG LAB and Fudan University.</p>
                        <el-divider></el-divider>
                        <p>If you used FACS, please cite:</p>
                        <p>Santos-Junior, C.D., Pan, S.J., Zhao, X.M., Coelho, L.P. 2019. FACS: antimicrobial peptide screening in genomes and metagenomes. BioRxiv xx.</p>
                        <el-divider></el-divider>
                        <h3>Contacts:</h3>
                        <p>Dr. Celio Dias Santos Junior (<a href="mailto:celio.diasjunior@gmail.com">celio.diasjunior@gmail.com</a>)</p>
                        <p>Dr. Luis Pedro Coelho (<a href="mailto:luis@luispedro.org">luis@luispedro.org</a>)</p>
                    </div>

                </div>
            </el-col>
            <el-col :span="2"><div class="grid-content"></div></el-col>
        </el-row>
    </div>
</template>

<script>
    export default {
        name: "About"
    }
</script>

<style scoped lang="less">
    .about-container{
        height: 100%;
        margin-top: 2em;
        /*background-color: greenyellow;*/
    }

    .content{
        border: 1px solid;
    }

    .content-header{
        border-bottom: 1px solid;
        padding: 1em 2em;
        background-color: antiquewhite;
    }

    .content-body{
        padding: 0.5em 2em;
    }

    p{
        font-family: Consolas;
        line-height: 1.5em;
        text-indent: 2em;
    }

</style>
