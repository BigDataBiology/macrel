<template>

    <el-container class="home-container">

        <el-header style="height: 3.6em">
            <div>
                <router-link :to="{name:'choice'}">
                    <el-button type="text"><h2 style="color: red">Macrel</h2></el-button>
                </router-link>
            </div>

            <div>
                <el-menu :default-active="$route.path" class="home-el-menu" mode="horizontal" @select="handleSelect" text-color="#ffffff" active-text-color="#20A0FF" background-color=transparent router >
                    <el-menu-item index="/prediction"><span>Prediction</span></el-menu-item>
                    <el-menu-item index="#"><a href="https://github.com/BigDataBiology/macrel" target="_blank">Download (Github)</a></el-menu-item>
                    <el-menu-item index="/about"><span>About</span></el-menu-item>
                    <el-menu-item index="/choice" hidden><span>choice</span></el-menu-item>
                </el-menu>
            </div>


            <div>
                <el-button style="padding: 0" circle><a href="https://github.com/BigDataBiology/macrel" target="_blank"><img src="../assets/image/Github_48px_1230449_easyicon.net.png"></a></el-button>
            </div>

        </el-header>

        <el-main style="margin-top: 0.1em;padding: 0.1em">
            <router-view/>
        </el-main>

        <el-footer style="height: 1.5rem"></el-footer>
    </el-container>

</template>

<script>
// @ is an alias to /src
// import HelloWorld from '@/components/HelloWorld.vue'

export default {
    name: 'home',

    data() {
        return {

        };
    },

    methods: {
        handleSelect(key, keyPath) {
            window.console.log(key, keyPath);
        },

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
        align-items: center;

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
