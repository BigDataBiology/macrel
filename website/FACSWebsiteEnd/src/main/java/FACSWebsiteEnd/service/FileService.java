package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.config.RemoteProperties;
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
    List readFileToObjectFromLocal(String fullFilePath, Object object);

    List readFileToObjectFromRemote(RemoteProperties remoteProperties,String filePath,Object object);
    FileInfo uploadFileToRemote(RemoteProperties remoteProperties,String dir,MultipartFile file);
    FileInfo saveContentToFileRemotely(RemoteProperties remoteProperties,String dir,String extension,String content);
}
