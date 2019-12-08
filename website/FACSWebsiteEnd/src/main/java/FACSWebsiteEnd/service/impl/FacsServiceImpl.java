package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FacsOutTsv;
import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.CommandUtils;
import FACSWebsiteEnd.utils.CommonUtils;
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
    public FileInfo saveSequenceToFile(String sequence, String dataType) {

        // 将序列文本输出为指定文件
        String extension = "";
        if (Constant.PEPTIDES.equals(dataType)){
            extension = Constant.FA;
        } else if (Constant.NUCLEOTIDE.equals(dataType)){
            extension = Constant.FASTTQ;
        }
        FileInfo fileInfo = fileService.saveTextToFile(sequence, extension);
        return fileInfo;

    }

    @Override
    public FileInfo saveFile(MultipartFile multipartFile) {
        FileInfo fileInfo = fileService.upload(multipartFile);
        return fileInfo;
    }

    @Override
    public List<Object> callShell(FileInfo fileInfo, String dataType) {

        String command = "";

        Map<String,Object> commandParams = new HashMap<String, Object>();
        String bash = Constant.BASH;
        String shellPath = Constant.FACS_HOME + Constant.FACS_SHELL;

        String outfolderName = CommonUtils.getCurrentTime();
        String outfolderPath = Constant.FACS_OUT_PARENT;

        if (Constant.PEPTIDES.equals(dataType)){
            commandParams.put("--mode","p");
            commandParams.put("--fasta",fileInfo.getFullpath());
            commandParams.put("-t",1);
            commandParams.put("--block",1000000);
            commandParams.put("--outfolder",outfolderPath);

            command = CommandUtils.buildShellCommand(bash,shellPath,commandParams);
//            System.out.println(command);

            // 远程执行
            //RemoteUtils.remoteInvokeShell(command);

            // 本地执行
            CommandUtils.executeLocalScript(command);

            // 读取tsv结果文件
            String fullFilePath = outfolderPath + "/" + Constant.FACS_OUT_FILENAME;
            List<Object> objects = fileService.readTsvGzToObject(fullFilePath, new FacsOutTsv());
            return objects;

        } else if(Constant.NUCLEOTIDE.equals(dataType)){
            // todo
            return null;
        }

        return null;
    }
}
