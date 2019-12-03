package FACSWebsiteEnd.utils;

import java.util.UUID;

/**
 * @Author: HiramHe
 * @Date: 2019/11/7 11:12
 */
public class CommonUtils {

    public static String getUUID(){
        return UUID.randomUUID().toString().replace("-","");
    }
}
