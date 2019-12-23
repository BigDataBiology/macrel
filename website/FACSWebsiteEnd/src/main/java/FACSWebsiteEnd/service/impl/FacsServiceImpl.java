package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.config.RemoteProperties;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommandUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:33
 * QQ:776748935
 */

@Service
public class FacsServiceImpl implements FacsService {

    @Autowired
    private FileService fileService;

    @Override
    public String callPipeline(String pipelineHome, FileInfo fileInfo, String currentOutDir, String dataType, RemoteProperties remoteProperties) {

        String command = "";

        Map<String,Object> commandParams = new HashMap<String, Object>();
        String bash = Constant.BASH;
        String shellPath = pipelineHome + Constant.FACS_SHELL;

        String tempFolderName = Constant.FACS_TEMP_FOLDER_PREFIX + fileInfo.getFilenameWithOutExtension();

        String inputFilePath = fileInfo.getPath();
        String outputFilePath = null;

        if (Constant.PEPTIDES.equals(dataType)){
            // run FACS on peptides

            commandParams.put("--mode","p");
            this.commonParams(commandParams,inputFilePath,currentOutDir,tempFolderName);

            command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);
//            System.out.println(command);

            this.execute(command,remoteProperties);

            outputFilePath = currentOutDir + Constant.FACS_OUT_FILENAME;

        } else if(Constant.CONTIGS.equals(dataType)){
            // run FACS on contigs

            commandParams.put("--mode","c");
            this.commonParams(commandParams,inputFilePath,currentOutDir,tempFolderName);

            command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);
//            System.out.println(command);

            this.execute(command,remoteProperties);

            outputFilePath = currentOutDir + Constant.FACS_OUT_FILENAME;

        }

        return outputFilePath;
    }

    private void commonParams(Map<String,Object> commandParams,String inputFilePath,String currentOutDir, String tempFolderName){
        commandParams.put("--fasta",inputFilePath);
        commandParams.put("-t",Constant.FACS_THREADS);
        commandParams.put("--block",Constant.FACS_BLOCK);
        commandParams.put("--outfolder",currentOutDir);
        commandParams.put("--tmp",tempFolderName);
    }

    private void execute(String command,RemoteProperties remoteProperties){
        if (!remoteProperties.getEnableRemote()){
            // 本地执行
            CommandUtils.executeCommandLocally(command);
        } else {
            // todo：远程执行
            CommandUtils.executeCommandRemotely(remoteProperties,command);
        }
    }

}
