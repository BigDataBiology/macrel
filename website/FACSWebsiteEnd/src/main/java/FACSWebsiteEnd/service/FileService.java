package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import org.springframework.web.multipart.MultipartFile;

import java.util.List;

/**
 * @Author: HiramHe
 * @Date: 2019/11/28 16:54
 * QQ:776748935
 */
public interface FileService {

    FileInfo uploadFileToLocal(MultipartFile file, String fullDir);

    FileInfo saveTextToFileLocally(String text, String fullDir, String extension);

    List<Object> readLocalTsvGzToObject(String fullFilePath, Object object);
}
