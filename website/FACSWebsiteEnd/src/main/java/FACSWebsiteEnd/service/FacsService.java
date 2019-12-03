package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.DataUploaded;
import FACSWebsiteEnd.Entity.FileInfo;
import org.springframework.web.multipart.MultipartFile;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 16:23
 * QQ:776748935
 */
public interface FacsService {

    FileInfo saveSequenceToFile(String sequence);
    FileInfo saveFile(MultipartFile multipartFile);
    Boolean callShell(String sequenceType, String mode, String read_1);
}
