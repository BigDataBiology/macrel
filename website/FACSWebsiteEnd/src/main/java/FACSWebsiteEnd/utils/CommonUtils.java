package FACSWebsiteEnd.utils;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.UUID;

/**
 * @Author: HiramHe
 * @Date: 2019/11/7 11:12
 */
public class CommonUtils {

    public static String getUUID(){
        return UUID.randomUUID().toString().replace("-","");
    }

    public static String getCurrentTime(){
        SimpleDateFormat sdf = new SimpleDateFormat();
        sdf.applyPattern("yyyy_MM_dd__HH_mm_ss");
        Date date = new Date();
        return sdf.format(date);
    }
}
