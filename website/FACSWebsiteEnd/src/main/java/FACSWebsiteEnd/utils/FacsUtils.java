package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.config.PipelineProperties;
import FACSWebsiteEnd.config.RemoteProperties;

/**
 * @Author: HiramHe
 * @Date: 2019/12/10 18:49
 * QQ:776748935
 */
public class FacsUtils {

    public static void createFolderLocally(String dir){
        // 在Linux上创建文件夹，用来保存用户上传的数据
        String space = " ";
        String commands = "mkdir -p" + space + dir;
        String[] shell = {"/bin/sh", "-c",commands};

        CommandUtils.executeCommandsLocally(shell);
    }

    public static String createFolderLocally(String dir, String folderName){
        // 在Linux上创建文件夹，用来保存当前结果
        String newDir = null;

        String space = " ";
        String commands = "cd" + space + dir + ";" +
                "mkdir -p" + space + folderName;
        String[] shell = {"/bin/sh", "-c",commands};

        CommandUtils.executeCommandsLocally(shell);

        newDir =  dir + folderName + "/";
        return newDir;
    }

    public static String createFolderRemotely(RemoteProperties remoteProperties,String dir,String folderName){

        String newDir = null;

        String space = " ";
        String commands = "mkdir -p" + space + dir + folderName + "/";

        CommandUtils.executeCommandRemotely(remoteProperties,commands);

        newDir = dir + folderName + "/";
        return newDir;
    }

    public static void createFolderRemotely(RemoteProperties remoteConfiguration, PipelineProperties pipelineConfiguration){

        String space = " ";
        String commands = "mkdir -p" + space + pipelineConfiguration.getInputDir() + ";" +
                "mkdir -p" + space + pipelineConfiguration.getOutputDir();

        CommandUtils.executeCommandRemotely(remoteConfiguration,commands);
    }
}
