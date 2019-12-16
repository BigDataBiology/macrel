package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import org.springframework.web.multipart.MultipartFile;

import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:23
 * QQ:776748935
 */
public interface FacsService {

    void callShellScript(FileInfo fileInfo, String currentOutDir, String dataType, Boolean isRemote);
    List<Object> readLocalResults(String filePath);
    List<Object> readRemoteResults(String outfolderPath,String filename);

}
