package FACSWebsiteEnd.service;

import FACSWebsiteEnd.Entity.FileInfo;
import FACSWebsiteEnd.common.ResultCode;
import FACSWebsiteEnd.common.ResultObject;
import FACSWebsiteEnd.config.RemoteProperties;
import org.springframework.web.multipart.MultipartFile;

import javax.servlet.http.HttpServletResponse;
import java.io.File;
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
    ResultObject getFileForDownloadFromLocal(File file, String filename, HttpServletResponse response);

    List readFileToObjectFromRemote(RemoteProperties remoteProperties,String filePath,Object object);
    FileInfo uploadFileToRemote(RemoteProperties remoteProperties,String dir,MultipartFile file);
    FileInfo saveContentToFileRemotely(RemoteProperties remoteProperties,String dir,String extension,String content);
    ResultObject getFileForDownloadFromRemote(RemoteProperties remoteProperties,String filePath,String filename,HttpServletResponse response);
}
