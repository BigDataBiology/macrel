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

    FileInfo saveSequenceToFile(String sequence, String dataType);
    FileInfo saveFile(MultipartFile multipartFile);
    List<Object> callShell(FileInfo fileInfo, String dataType);

}
