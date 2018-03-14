import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

export class WelcomePage extends React.Component<{},{}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    /**
     * Render the component.
     */
    render() {


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
                            
                            <ul>
                                <li>The database currently consists of organisms</li>
                                <li>The database currently consists of organisms</li>
                                <li>The database currently consists of organisms</li>
                            </ul>
        
                        </CardText>
                    </Card>
            </div>

        );
    }

}