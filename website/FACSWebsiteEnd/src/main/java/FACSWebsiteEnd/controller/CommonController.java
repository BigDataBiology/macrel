package FACSWebsiteEnd.controller;

import FACSWebsiteEnd.common.ResultObject;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

/**
 * @Author: HiramHe
 * @Date: 2019/12/8 15:39
 * QQ:776748935
 */

@RestController
@CrossOrigin
@RequestMapping("/common")
public class CommonController {

    @GetMapping("/testCommon")
    public ResultObject testCommon(){
        return ResultObject.success();
    }
}
