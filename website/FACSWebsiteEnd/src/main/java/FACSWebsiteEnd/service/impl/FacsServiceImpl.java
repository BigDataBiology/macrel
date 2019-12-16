package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommandUtils;
import FACSWebsiteEnd.utils.CommonUtils;
import FACSWebsiteEnd.utils.FacsUtils;
import FACSWebsiteEnd.utils.RemoteUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.web.multipart.MultipartFile;

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
    FileService fileService;

    @Override
    public void callShellScript(FileInfo fileInfo, String currentOutDir, String dataType, Boolean isRemote) {

        String command = "";

        Map<String,Object> commandParams = new HashMap<String, Object>();
        String bash = Constant.BASH;
        String shellPath = Constant.FACS_HOME + Constant.FACS_SHELL;

        String tempFolderName = Constant.FACS_TEMP_FOLDER_PREFIX + fileInfo.getFilenameWithOutExtension();

        if (Constant.PEPTIDES.equals(dataType)){
            // run FACS on peptides

            commandParams.put("--mode","p");
            commandParams.put("--fasta",fileInfo.getPath());
            commandParams.put("-t",1);
            commandParams.put("--block",1000000);
            commandParams.put("--outfolder",currentOutDir);
            commandParams.put("--tmp",tempFolderName);

            command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);
//            System.out.println(command);

            if (!isRemote){
                // 本地执行
                CommandUtils.executeLocalCommand(command);
            }else {
                // todo：远程执行
                RemoteUtils.remoteInvokeShell(command);
            }

        } else if(Constant.CONTIGS.equals(dataType)){
            // run FACS on contigs

            commandParams.put("--mode","c");
            commandParams.put("--fasta",fileInfo.getPath());
            commandParams.put("-t",1);
            commandParams.put("--block",100000);
            commandParams.put("--outfolder",currentOutDir);
            commandParams.put("--tmp",tempFolderName);

            command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);
//            System.out.println(command);

            if (!isRemote){
                // 本地执行
                CommandUtils.executeLocalCommand(command);
            } else {
                // todo：远程执行
                RemoteUtils.remoteInvokeShell(command);
            }

        }
    }

    @Override
    public List<Object> readLocalResults(String filePath) {

        // 读取tsv结果文件
        List<Object> objects = fileService.readLocalTsvGzToObject(filePath, new FacsOutTsv());
        return objects;
    }

    @Override
    public List<Object> readRemoteResults(String outfolderPath, String filename) {
        // todo
        return null;
    }
}
