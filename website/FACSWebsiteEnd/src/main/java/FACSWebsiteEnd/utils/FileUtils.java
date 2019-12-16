package FACSWebsiteEnd.utils;

import org.springframework.web.multipart.MultipartFile;

import java.util.HashMap;
import java.util.Map;

/**
 * @Author: HiramHe
 * @Date: 2019/12/2 17:11
 * QQ:776748935
 */
public class FileUtils {

    public static Map getFileInformation(MultipartFile file){

        Map<String,String> infoMap = new HashMap<String,String>();

        String originalFilename = file.getOriginalFilename();

        // 获取上传文件的类型
        String extension = originalFilename.indexOf(".") != -1 ?
                originalFilename.substring(originalFilename.lastIndexOf(".")+1) : null;

        int lengthOfSuffix = extension.length() + 1;

        String filenameWithOutExtension = originalFilename.substring(0,originalFilename.length()-lengthOfSuffix);

        infoMap.put("filenameWithExtension",originalFilename);
        infoMap.put("filenameWithOutExtension",filenameWithOutExtension);
        infoMap.put("extension",extension);

        return infoMap;
    }
}
