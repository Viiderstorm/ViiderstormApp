require.config(
    {
        paths: {
            text: 'text',
            knockout: 'knockout',
            jquery: 'jquery',
            three: '../three',
            semantic: '../semantic'
        },
        shim: {
            'semantic': ['jquery'],
            'three':{
                exports: 'THREE'
            }
        }
    }
)


define('main', function(){
    
        require(['knockout', 'jquery', 'semantic', 'text'], function(ko, $, semantic){
            $(".launch.button").click(function(){
                    $(".ui.sidebar").sidebar('toggle')
            });
            console.log("Done Loading");
        })
        
    }
)