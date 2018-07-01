

export default class config {

    static restServer = 'http://localhost'
    static restPort = '5001'

    static getRestAddress()
    {
        return this.restServer + ":" + this.restPort ;
    }

    static table_even_bg = 'lightgrey';

    static axiosConfig = {
        crossdomain: true,
        headers: {
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET, PUT, POST, DELETE, OPTIONS"
        }
      };
 
}