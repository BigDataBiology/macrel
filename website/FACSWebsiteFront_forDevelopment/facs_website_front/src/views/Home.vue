<template>

    <el-container class="home-container">
        <!--页面头部区-->
        <el-header style="height: 3.6em">
            <div><h2>FACS</h2></div>

            <div>
                <el-menu :default-active="activeIndex" class="home-el-menu" mode="horizontal" @select="handleSelect" text-color="#ffffff" active-text-color="#20A0FF" background-color=transparent router >
                    <el-menu-item index="/prediction"><span>Prediction</span></el-menu-item>
                    <el-menu-item index="#"><a href="https://github.com/FACS-Antimicrobial-Peptides-Prospection/FACS" target="_blank">Github</a></el-menu-item>
                    <el-menu-item index="/about"><span>About</span></el-menu-item>
                </el-menu>
            </div>

            <div>
                <el-button style="padding: 0" circle><a href="https://github.com/FACS-Antimicrobial-Peptides-Prospection/FACS" target="_blank"><img src="../assets/image/Github_48px_1230449_easyicon.net.png"></a></el-button>
            </div>

        </el-header>

        <!--主体区-->
        <el-main style="margin-top: 0.1em;padding: 0.1em">
            <router-view/>
        </el-main>

        <!--尾部区-->
        <el-footer style="height: 1.5rem"></el-footer>
    </el-container>

</template>

<script>
// @ is an alias to /src
// import HelloWorld from '@/components/HelloWorld.vue'

export default {
    name: 'home',

    methods: {
        handleSelect(key, keyPath) {
            console.log(key, keyPath);
        }
    },

    data() {
        return {
            activeIndex: "/prediction",
        };
    },

}
</script>

<style lang="less" scoped>

    .el-container{
        height: 100%;
    }

    .el-header{

        color: red;
        background-color: #34495E;
        display: flex;
        justify-content: space-evenly;
        /* 内容垂直居中 */
        align-items: center;

        /* 嵌套布局 */
        > div {
            display: flex;
            align-items: center;

            .el-menu{
                border-bottom: none;
                .el-menu-item{
                    font-size: 1.3em;
                    padding: 0 2em;
                }

            }
        }
    }

    a{
        text-decoration: none;
    }

    .el-footer{
        height: 30rem;
    }


</style>
