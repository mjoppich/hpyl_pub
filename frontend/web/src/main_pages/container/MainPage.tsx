import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';
import axios from 'axios';
import config from '../config';

export class WelcomePage extends React.Component<{},{stats: any}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentWillMount() {
        var self=this;
        self.setState({stats: {}})

        axios.get(config.getRestAddress() + "/stats", config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var stats = response.data;
            self.setState({stats: stats})

        })
        .catch(function (error) {
            console.log(error)
            self.setState({stats: {}})
        });
    }

    /**
     * Render the component.
     */
    render() {

        var statElems = <p>No database stats available</p>;

        if (this.state.stats.org_count)
        {
            statElems = <ul>
            <li>The database currently consists of {this.state.stats.org_count} organisms.</li>
            <li>The database currently includes {this.state.stats.hom_count} simple homologies.</li>
            <li>The database currently includes {this.state.stats.comb_count} combined homologies.</li>
            <li>The database currently includes {this.state.stats.mul_comb_count} multiple-combined homologies.</li>
            </ul>;
        }
        return (

            <div>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="Welcome!"
                    subtitle="heliPyloriBD v0.1"
                />
                <CardText >
                    <p>heliPyloriDB v.01 online.</p>
                    <p>Explore the database in the explore tab</p>
                    <p>Provide some statistics of the project here</p>

                </CardText>
            </Card>
                        <Card>
                        <CardTitle
                            title="heliPyloriBD v0.1"
                            subtitle="Statistics"
                        />
                        <CardText >
                            
                            {statElems}
        
                        </CardText>
                    </Card>
            </div>

        );
    }

}