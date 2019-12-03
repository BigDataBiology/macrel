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

        String type = null;
        String filename = file.getOriginalFilename();

        infoMap.put("filename",filename);

        // 获取上传文件的类型
        type = filename.indexOf(".") != -1 ?
                filename.substring(filename.lastIndexOf(".")+1) : null;
        infoMap.put("type",type);

        return infoMap;
    }
}
