package FACSWebsiteEnd.utils;

import FACSWebsiteEnd.common.Constant;

/**
 * @Author: HiramHe
 * @Date: 2019/12/10 18:49
 * QQ:776748935
 */
public class FacsUtils {

    public static String makeSavedFolderOnLinux(String jarPath){
        // 在Linux上创建文件夹，用来保存用户上传的数据
        String commands = "cd "+jarPath+";"+"mkdir "+ Constant.FILESAVED_FOLDER_NAME;
        String[] shell = {"/bin/sh", "-c",commands};
        CommandUtils.executeLocalCommandArray(shell);

        String savedDir = jarPath + "/" + Constant.FILESAVED_FOLDER_NAME + "/";
        return savedDir;
    }

    public static String makeAllOutFolderOnLinux(String jarPath){
        // 在Linux上创建文件夹，用来保存所有结果
        String commands = "cd "+jarPath+";"+"mkdir "+Constant.FACS_ALLOUT_FOLDER_NAME;
        String[] shell = {"/bin/sh", "-c",commands};
        CommandUtils.executeLocalCommandArray(shell);

        String allOutDir = jarPath +  "/" + Constant.FACS_ALLOUT_FOLDER_NAME + "/";
        return allOutDir;
    }

    public static String makeCurrentOutFolderOnLinux(String allOutDir, String foldername){

        // 在Linux上创建文件夹，用来保存当前结果
        String commands = "cd "+allOutDir+";"+"mkdir "+foldername;
        String[] shell = {"/bin/sh", "-c",commands};
        CommandUtils.executeLocalCommandArray(shell);

        String currentOutDir =  allOutDir + foldername + "/";
        return currentOutDir;

    }
}
