package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import org.springframework.web.multipart.MultipartFile;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:54
 * QQ:776748935
 */
public interface FileService {

    FileInfo upload(MultipartFile file);

    FileInfo saveTextToFile(String text, String extension);
}
