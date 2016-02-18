require.config(
    {
        paths: {
            sailsio: '/js/dependencies/sails.io',
            jquery: '/js/jquery',
            three: '/js/three',
            semantic: '/js/semantic'
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
    
        require(['sailsio', 'jquery', 'semantic'], function(sails, $, semantic){
            console.log("Done Loading");  
        })
        
    }
)