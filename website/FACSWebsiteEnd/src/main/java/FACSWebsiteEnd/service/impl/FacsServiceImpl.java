package FACSWebsiteEnd.service.impl;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.Constant;
import FACSWebsiteEnd.service.FacsService;
import FACSWebsiteEnd.service.FileService;
import FACSWebsiteEnd.utils.RemoteUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.util.StringUtils;
import org.springframework.web.multipart.MultipartFile;

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
    public FileInfo saveSequenceToFile(String sequence) {

        // 将序列文本输出为指定文件
        String extension = Constant.EXTENSION;
        FileInfo fileInfo = fileService.saveTextToFile(sequence, extension);
        return fileInfo;

    }

    @Override
    public FileInfo saveFile(MultipartFile multipartFile) {
        FileInfo fileInfo = fileService.upload(multipartFile);
        return fileInfo;
    }

    @Override
    public Boolean callShell(String sequenceType, String mode, String read_1) {

        String ip = "39.106.68.204";
        int port = 22;
        String username = "HiramHe";
        String password = "hiram1024";

        String space = " ";

        String command1 = "bash"
                +space+"helloWorld04.sh"
                +space+sequenceType
                +space+mode
                +space+read_1;

        RemoteUtils.remoteInvokeShell(ip,port,username,password,command1);

        return null;
    }
}
