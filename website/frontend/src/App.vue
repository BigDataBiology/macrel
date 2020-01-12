<template>
    <div id="app">
        <!-- placeholder for router -->
        <router-view />
    </div>
</template>

<script>

export default {
    name: 'app'
}
</script>

<style scoped lang="less">

    html{
        height:100%;
    }

    body{
        min-height:100%;
        margin:0;
        padding:0;
        position:relative;
    }

    .el-container{
      /*height: 100vh;*/
    }
    .el-header {
        /*background-color: #B3C0D1;*/
        /*color: #333;*/
        text-align: center;
        line-height: 60px;
    }

  .el-main {
    /*background-color: #E9EEF3;*/
    color: #333;
    /*text-align: center;*/
    /*line-height: 160px;*/
  }

    .el-footer{
        text-align: center;
        line-height: 60px;
        position: absolute;
        bottom:0;
        width:100%;height:100px;
    }

    .el-col {
        border-radius: 4px;
    }

    .grid-content {
        border-radius: 4px;
        min-height: 36px;
    }

    .router-link-active{
        text-decoration: none;
    }

</style>
